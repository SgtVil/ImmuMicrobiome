#' Plot the overall variance of your dataset
#'
#' @description
#' This function will plot the overall variance explained by the selected factors.
#' You can choose to plot it as boxplots or heatmap.
#' - Boxplots : points will represent the mean abundance of the features. Features with p.adj < 0.05 and p < 0.05 will be plotted.
#' - Heatmap : all factors will be plotted as rows. Currently doesn't shows which features ar significant.
#'
#' @param variance An object returned by the \link{calculate_variance}
#' @param top The top N features to plot. Default = 30
#' @param plot_type Categorical. Plot as boxplots or heatmap. Default= "boxplot"
#'
#' @return
#' @export
#'
#' @examples
plot_all_variance = function(variance, top=30, plot_type= "boxplot",  col= c("brown", "orange", "grey")){
  var.exp = variance$variance %>%
    filter(variable=="var.exp")%>%
    pivot_longer(names_to = "factor", values_to = "value", cols = !features:variable)%>%
    full_join(variance$p.value, by=c("features","factor"))%>%
    full_join(variance$variance %>%
                filter(variable=="mean.feat") %>%
                pivot_longer(names_to = "factor", values_to = "mean.feat", cols = !features:variable),
              by=c("features", "factor"))
  if(plot_type=="boxplot"){

    p =var.exp %>%
      ggplot(aes(y = value, fct_reorder(factor, var.exp%>%
                                          group_by(factor)%>%
                                          mutate(mean.fac = median(value, na.rm=T))%>%
                                          pull(mean.fac))))+

      geom_jitter(aes(size=mean.feat, fill= ifelse(p.adj<0.05, col[1],
                                                   ifelse(p<0.05, col[2], col[3]))),
                  shape=21)+
      geom_boxplot(size=1.5, alpha=0)+
      scale_fill_identity(name = "P values", breaks = c(col),
                          labels = c("FDR < 0.05", "nFDR < 0.05", "NS"),
                          guide = "legend")+
      guides(size= guide_legend("Mean features"))+
      coord_flip()+
      theme_bw()+
      labs(y="Variance explained by factors")+
      theme(axis.text.y = element_text(face="bold"),
            axis.title.y = element_blank(),
            legend.position = "bottom")
  }
  if(plot_type=="heatmap"){

    top= variance$variance %>%
      filter(variable=="var.tot")%>%
      pivot_longer(names_to = "factor", values_to = "value", cols = !features:variable)%>%
      group_by(features)%>%
      summarise(mean.fac= mean(value, na.rm=T)) %>%
      top_n(50, wt = mean.fac)

    p = var.exp %>%
      filter(features %in% top$features)%>%
      ggplot(aes(fct_reorder(features, mean.feat), factor))+
      geom_tile(aes(fill=value))+
      scale_fill_gradient(high = col[1], low = "white", name="Percent of variance")+
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
            axis.title = element_blank(),
            legend.position = "top")

  }
return(p)

}

