#' Plot volcano
#'
#' Make a volcano plot
#' @param volcano An object returned by the \link{fold_change} function
#' @param col Color vector.
#' @param ... Addtionnal arguments passed to \link{geom_point(aes(...))}
#'
#' @return
#' @export
#'
#' @examples
plot_volcano = function(volcano, size=2, col= c("#771155", "#77CCCC","grey"), ...){


  labels = volcano%>%
    filter(p<0.05)
  max.fold = max(abs(volcano$fold))
  p= volcano %>%
    mutate(col.val =ifelse( volcano %>%
                              pull(p.adj)<0.05, col[1], ifelse( volcano %>%
                                                                  pull(p)<0.05, col[2], col[3])))%>%
    ggplot(aes(fold, -log10(p.adj)))+
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
