#' Make circular plots for variance analysis.
#'
#' @description
#' A ggplot2 based function to plot variance explained by factors in a circular manner.
#'
#' @param variance Object returned by the function \link{calculate_variance}.
#' @param adjust Value to adjust the labels
#' @param show.percentage Logical. Show percentage for each features.
#' @param label.size Size of the labels.
#' @param text.size Size of the labels
#'
#' @return
#' A ggplot2 graphic.
#' @export
#'
#' @examples
#'
#' var = calculate_variance(metabolomic[,1:4], features = metabolomic[,5:217], cores = 5)
#' circular_variance_plot(variance, adjust=0.4, show.percentage = T)
circular_variance_plot= function(variance, adjust=1, col= c("brown","darkgreen","orange","violet"), show.percentage=F, label.size=5, text.size= 5){
  var.exp = variance$variance %>%
    filter(variable=="var.exp")%>%
    pivot_longer(names_to = "factor", values_to = "value", cols = !features:variable)%>%
    full_join(variance$p.value, by=c("features","factor"))%>%
    full_join(variance$variance %>%
                filter(variable=="mean.feat") %>%
                pivot_longer(names_to = "factor", values_to = "mean.feat", cols = !features:variable),
              by=c("features", "factor")) %>%
    mutate(features= sort(features))



  tmp=var.exp %>%
    # group_by(features)%>%
    mutate(features= sort(features))%>%
    summarise(
      id = 1:length(unique(var.exp$features)),
      angle = 90 - 360 * (1:length(unique(var.exp$features))-0.5) /
        length( unique(var.exp$features)),
      angle = ifelse(angle < -90, angle + 180, angle),
      label = unique(var.exp$features),
      hjust = ifelse(angle < -90, 1, 0),
      value = var.exp %>%
        group_by(features)%>%
        mutate(features= sort(features),
               value= sum(value)) %>%
        distinct(value)%>%
        pull(value)
    )

  p= var.exp %>%
    arrange(features)%>%
    # mutate()
    ggplot( aes(features, value))+
    geom_bar(aes( fill=factor),stat = "identity", color="black")+
    scale_fill_manual(values=col)+
    guides(fill=guide_legend("Explained\nvariance"))+
    geom_text(data= tmp, aes(label=label, as.numeric(factor(label)),  angle=angle, y= value+adjust), size=text.size)+
    coord_polar(start=0)+
    theme_minimal()+
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.text.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.margin = unit(c(-0.5, -0.5, 0, -0.5), "cm")
    )
  if(isTRUE(show.percentage)){
    p= p +
      geom_label(data= var.exp, aes(label= paste(round(value*100, 1), "%")), position =position_stack(0.5),
      size= 3, show.legend = F)
  }

 return(p)
}

