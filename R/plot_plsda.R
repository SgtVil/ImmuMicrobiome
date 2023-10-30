#' Plot PLS-DA using mixOmics `plsda`
#'
#' @description
#' A ggplot2 based function to plot PLS-DA based on the object returned by \code{\link[mixOmics]{plsda}}.
#' The function will also plot the features loadings, a biplot then. This function can also handle results of a `splsda`.
#'
#' @param res.plsda An object returned by \code{\link[mixOmics]{plsda}}.
#' @param comp Which component to plot.
#' @param color_vector Color vector
#' @param shape Shape of the dots.
#' @param top.loads Top N features to plot. Default = 10.
#' Loads are defined by the absolute value of each features for each components. If `top.loads = 10` the function
#' will keep the top10 features of component 1 and the top 10 features of the component 2.
#'
#' @return
#'  A `ggplot2` object.
#' @export
#'
#' @examples
#'
#' data(metabolomic)
#' res= mixOmics::plsda(X = metabolomic[,5:217], Y = metabolomic$birth_type, ncomp=5)
#'
#' # Basic plotting
#' plot_plsda(res, comp= c(1,2), shape=21, top.loads= 20)
#'
#' # Improve a bit the plot.
#' plot_plsda(res, comp= c(1,2), shape=21, top.loads= 20, color_vector= tol21rainbow[c(2,14)])+
#' guides(fill= guide_legend("Birth route"))+
#' theme_bw()+
#' theme(legend.position = "bottom",
#'        axis.title = element_text(size=15, face="bold"),
#'        panel.grid = element_blank(),
#'        legend.text = element_text(size=15),
#'        legend.title = element_text(size=20, face="bold"),
#'        axis.text = element_text(size=15))
#'
#'


plot_plsda = function(res.plsda, comp= c(1, 2), color_vector=  tol21rainbow, shape, top.loads= 10){
  loads = res.plsda$loadings$X%>%
    as.data.frame()%>%
    arrange(desc(abs(comp1))) %>%
    top_n(top.loads, wt= abs(comp1)) %>%
    rownames_to_column()%>%
    full_join(res.plsda$loadings$X%>%
                as.data.frame()%>%
                arrange(abs(desc(comp2))) %>%
                top_n(top.loads, wt= abs(comp2))%>%
                rownames_to_column())
  scaler <- max(res.plsda$variates$X, na.rm = TRUE)/ max(abs(res.plsda$loadings$X), na.rm = TRUE)
  variates = res.plsda$variates$X %>%
    bind_cols(res.plsda$Y)%>%
    as.data.frame()

  var.exp = round(res.plsda$prop_expl_var$X*100, 1)

  variates %>%
  ggplot(aes(x = variates[, comp[1]],
              y = variates[, comp[2]]))+
    geom_point(shape=21, aes(fill=variates[,length(colnames(variates))]), size=3)+
    stat_ellipse(aes(fill=variates[,length(colnames(variates))]), color="black", geom = "polygon", alpha=0.3)+
    geom_segment(data = loads , aes(xend=comp1*scaler*2, yend=comp2*scaler*2, x=0, y=0), alpha=0.3,
                 arrow = arrow(angle = 30, length  = unit(0.25, "cm"),
                               ends = "last", type = "open"))+
    geom_label_repel(data= loads, aes(label=rowname, x=comp1*scaler*2, y=comp2*scaler*2))+
    scale_fill_manual(values = color_vector)+
    labs(x=paste("Component 1:", var.exp[1],"%"),
         y= paste("Component 2:", var.exp[2], "%"))
}
