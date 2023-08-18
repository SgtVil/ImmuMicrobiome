plot_bubble_ig_seq = function(igseq, taxa_rank="Genus", filtering_taxa=5, grouping_factor){


  tmp= igseq %>%
    group_by(Genus, mother_infant, timepoint)%>%
    mutate(n=n())%>%
    filter(n>filtering_taxa , Genus!="NA")%>%
    wilcox_test(score~1, mu=0, alternative ="two.sided", detailed = T)

  tmp2= igseq %>%
    group_by(Phylum, Genus, mother_infant, timepoint)%>%
    mutate(n=n())%>%
    filter(n>filtering_taxa, Genus!="NA")%>%
    dplyr::select(Genus, score, unique_timepoint)%>%
    mutate(mean_score= median(score))%>%
    left_join(tmp)

  tmp2$mean_score[tmp2$mean_score <= (-5)]= (-5)
  tmp2$unique_timepoint= factor(tmp2$unique_timepoint, levels = c("child_delivery","child_2_months","child_24_months", "mother_24_months"))

  tmp2= tmp2%>%
    ungroup()%>%
    select(p, score, mean_score, Phylum, Genus, mother_infant, unique_timepoint)%>%
    mutate(alpha= ifelse(tmp2$p<0.05, -log10(tmp2$p), 0),
           alpha2= ifelse(tmp2$mean_score<0, -alpha, alpha),
           mean_score2= rescale(c(abs(mean_score)), to=c(0,5)))

  p1= tmp2 %>%
    ggplot(aes(0, Genus,  fill=alpha2, size=mean_score2)) +
    geom_point(shape=21) +
    scale_fill_gradient2(low = "olivedrab4", mid = "white", high = "brown", midpoint = 0, guide=F)+
    scale_size_continuous(breaks=c(0,3,3,4,4,5,5),labels=c(-5,-4,-3,0,3,4,5),range = c(0,5))+
    guides(size = guide_legend( override.aes = list(fill =colorRampPalette(c("darkgreen", "white", "brown"))(7),
                                                    size=c(5,4,3,1,3,4,5)), nrow=1,
                                direction = "horizontal",
                                title.position = "top",
                                label.position = "bottom",
                                label.hjust = 0.5,
                                label.vjust = 1))+
    labs(fill="", x="Age", size="IgAseq median score")+
    facet_grid(Phylum~unique_timepoint, scales="free", space = "free", labeller = labeller(unique_timepoint=fac))+
    theme(axis.text.y = element_text(size=8, face="bold"),
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank(),
          strip.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = 'bottom',
          legend.title = element_text(size=10, face="bold"),
          legend.text = element_text(size=10, face="bold"),
          legend.box = "vertical",
          legend.box.background =  element_rect(colour = "black"),
          strip.background = element_blank(),
          strip.text = element_text(size=12, face="bold"),
          panel.grid.major.y = element_line(linetype = 2, size=.05)
    )
}
