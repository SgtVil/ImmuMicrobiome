plot_ig_seq = function(ig_seq, type= c("plotly", "ggplot"), save=F){
  
    df= ig_seq@ig_seq_all
    df = mutate(df, tooltip = paste(taxonomy, "\nslide_z: ",
                                          round(score, digits=3)))
    ell= ig_seq@ellipse_data
    ell= ell %>% mutate(tooltip = df$tooltip)
    stat = df$score>=1.96 | df$score <= -1.96
  
    p <- ggplot(df)+
      geom_polygon(data = ell, aes(ell[,1], ell[,2]), alpha=0.2, fill="burlywood", linetype=2)+
      geom_path(data = ell, aes(ell[,3], ell[,4]), color="chocolate1", linetype=2)+
      geom_path(data = ell, aes(ell[,5], ell[,6]), color="coral3", linetype=2)+
      geom_point(aes(log10_abundance, 
                     log2_ratio,
                     size=ellipse_level, 
                     color=as.character(ellipse_level), 
                     text=tooltip, 
                     shape=stat),                       alpha=0.6) +
      geom_point(aes(log10_neg_abundance, log2_neg_ratio), color="black")+
      scale_color_manual(values=c("cyan4","burlywood","chocolate1", "coral3"), name= "Confidence interval")+
      labs(x= "Log10 pos * neg", y="Log2 pos/neg")+
      theme_minimal()
    
    if(type=="ggplot"){
      if(save==T){
        ggsave(filename = paste0("./IgAseq plots/png/", unique(df[,"sample_id"]),
                                           "_",zero_treatment ,'_.png'), plot = p)
        dir.create("./IgAseq plots/png", recursive = T)
        
        return(p)
      } else {
        
        return(p)
        }
     
    }
    
    if(type=="plotly"){
      pp = plotly::ggplotly(p, tooltip ="text")
      if(save==T){
        
        saveWidget(pp, file=paste0("./IgAseq plots/plotly/", 
                                   unique(df[,"sample_id"]),
                                   colnames(df)[2], "_vs_", 
                                   colnames(df)[3], "_vs_", 
                                   colnames(df)[4],"_", 
                                   zero_treatment, 
                                   "_ratio_intensity_ellipse.html"))
        dir.create("./IgAseq plots/plotly", recursive = T)
        return(pp)
      } else {
       
        return(pp)
      }
     
    }
    
    
  
}
