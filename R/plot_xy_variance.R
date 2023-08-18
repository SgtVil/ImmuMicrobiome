#' Plot the variance for two factors.
#'
#' @param variance An object returned by \link{calculate_variance}.
#' @param x Factor to be plotted as x.
#' @param y Factor to be plotted as y.
#' @param fac The factor to depict p.values. Only the points with p.adj < 0.05
#' and p < 0.05 will be colorized and labelled
#' @param col Color vector.
#'
#' @return
#' @export
#'
#' @examples
plot_xy_variance = function(variance, x, y, fac, col= c("brown", "orange", "grey"), max.overlap=20){



  labels = variance$p.value %>%
    filter(factor==fac & p<0.05)

  max.var= variance$variance %>%
    filter(variable=="var.exp")%>%
    summarise(across(c(!!sym(x), !!sym(y)), .fns=max))%>%
    max
  min.var= variance$variance %>%
    filter(variable=="var.exp")%>%
    summarise(across(c(!!sym(x), !!sym(y)), .fns=min))%>%
    min
  p=variance$variance %>%
    filter(variable=="var.exp")%>%
    mutate(col.val= ifelse( variance$p.value %>%
                     filter(factor==fac)%>%
                     pull(p.adj)<0.05, col[1], ifelse( variance$p.value %>%
                                                         filter(factor==fac)%>%
                                                         pull(p)<0.05, col[2], col[3])))%>%
    ggplot()+
    geom_point(aes(fill=col.val,!!sym(x), !!sym(y), size=variance$mean.feat$mean.feat), shape=21)+
    scale_fill_identity(name = "P values", breaks = c(col),
                        labels = c("FDR < 0.05", "nFDR < 0.05", "NS"),
                        guide = "legend")+
    geom_abline(slope = 1, intercept = 0, linetype="dashed")+
    geom_label_repel(data= variance$variance %>%
                       filter(features %in% labels$features & variable=="var.exp"),
                     aes(label= features, !!sym(x), !!sym(y)), max.overlaps = max.overlap)+
    guides(size= guide_legend("Mean features"))+
    scale_y_continuous(limits = c(min.var, max.var), expand=c(0,0.005), labels=scales::percent)+
    scale_x_continuous(limits = c(min.var, max.var), expand=c(0,0.005), labels=scales::percent)
return(p)

}

