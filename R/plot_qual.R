#' Check the length of your reads.
#'
#' @param list An object returned by \link{summarise_fastq}
#' @param color_vector A color vector
#' @param n Number of fastq you used in \link{summarise_fastq}
#' @return
#' A ggplot.
#' @export
#' @import ggridges
#' @examples
#' TBD
plot_qual= function(list, n=5, color_vector=NULL, project_names){
  if(is.null(color_vector)) color_vector= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD",
                                            "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77",
                                            "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411",
                                            "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")


  df = purrr::list_transpose(list)

  depth= df$depth %>%
    transpose
p1=  do.call("rbind", depth$fwd)%>%
    as.data.frame()%>%
    rownames_to_column()%>%
    pivot_longer(cols=2:all_of(n)) %>%
    ggplot(aes(rowname, value, fill=rowname))+
    geom_boxplot()+
  labs(y="Sequencing depth")+
  scale_fill_manual(values=color_vector)+
  scale_color_manual(values= color_vector)+
  scale_y_log10()+
theme(axis.text.x =  element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank())

quality= df$quality %>%
  transpose
p2=  do.call("rbind", quality$fwd)%>%
  as.data.frame()%>%
  rownames_to_column()%>%
    mutate(project=  str_split(string = rowname, pattern = "\\.", simplify = T)[,1]) %>%
  ggplot(aes( mean, y=project, color=project, fill=project))+
  ggridges::geom_density_ridges2()+
  labs(x="Mean quality")+
  scale_fill_manual(values=color_vector)+
  scale_color_manual(values= color_vector)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

length= df$length %>%
  transpose
p3= do.call("rbind", length$fwd) %>%
  as.data.frame()%>%
  rownames_to_column()%>%
  pivot_longer(cols=2:1001)%>%
  ggplot(aes(x=rowname, value, fill=rowname, color=rowname))+
  geom_boxplot()+
  geom_jitter()+
  # geom_density()
  # ggridges::geom_density_ridges2()+
  labs(y="Read length")+
  scale_fill_manual(values=color_vector)+
  scale_color_manual(values= color_vector)+
  theme(axis.text.x =  element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
p= ggpubr::ggarrange(p1 ,p2, p3, common.legend = T, ncol=3, align="hv")
return(p)
}
