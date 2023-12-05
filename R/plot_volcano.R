#' Plot volcano
#'
#' Make a volcano plot based on the object returned by \link{fold_change}
#' @param volcano An object returned by the \link{fold_change} function
#' @param col Color vector.
#' @param ... Addtionnal arguments passed to \link{geom_point(aes(...))}
#'
#' @return
#' A ggplot2 object with features labelled depending on the p-values and the corrected p-values.
#' @export
#'
#' @examples
#'fold= fold_change(clinical_data = metabolomic[,1:4], metabolomic[,6:217], cores = 10)
#' lapply(fold, plot_volcano)
#' plot_volcano(fold$MeCT, size=2)

plot_volcano = function(volcano, size=2, col= c("#771155", "#77CCCC","grey"), ...){


  labels = volcano%>%
    filter(p<0.05)
  max.fold = max(abs(volcano$fold))
  p= volcano %>%
    mutate(col.val =ifelse( volcano %>%
                              pull(p.adj)<0.05, col[1], ifelse( volcano %>%
                                                                  pull(p)<0.05, col[2], col[3])))%>%
    ggplot(aes(fold, -log10(p)))+
    geom_point(aes(fill= col.val, ...), size=size, shape=21)+
    geom_hline(yintercept = 1.3)+
    scale_fill_identity(name = "P values",
                        breaks = c(col),
                        labels = c("FDR < 0.05", "nFDR < 0.05", "NS"),
                        guide = "legend")+
      geom_label_repel(data = volcano%>%
                         filter(features %in% labels$features)
                       , aes(label=features))+
    scale_x_continuous(limits=c(-max.fold, max.fold))+
    labs(x=paste("Fold change", unique(volcano$comp)), )
  return(p)
}
