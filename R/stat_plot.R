stat_plot <- function(data, condition, value, title=NULL, labels= NULL,
                      geom= "boxplot", test="wilcox.test", signif.size= 5, hide.ns=T, label="p.format", paired=F){
  require(ggplot2)
  require(ggpubr)
  require(mgsub)
  old <- data[,condition]
  theme_set(theme_bw()+theme(axis.text = element_text(size = 15), axis.title = element_text(size=20),
                             panel.grid = element_blank()))

  if(anyNA(unique(data[,condition]))){
      data[,condition][is.na(data[,condition])] <- "NA"
        cat(unique(data[,condition]), "\n")}

  if(is.numeric(data[,condition])){
    data[,condition]<- as.factor(data[,condition])#get a hold on numeric data such as 0 and 1
  }
  map= aes_string(x=condition, y=value, fill=condition)#create aesthetics
  comp= combn(as.character(na.omit(unique(data[,condition]))), 2, simplify = F)#create combinations but without NAs

    if(anyNA(unique(data[,condition]))){
      warning("NA detected in your conditions, NA won't be computed")#warnings for NA presence
    }
  plot_type = switch(geom, "boxplot" = geom_boxplot(alpha=0.2),
                     "violin"= geom_violin(alpha=0.2),
                     "col"= geom_histogram(alpha=0.2, stat = "identity"),
                     'jitter'= geom_blank(),
                     "crossbar"=geom_crossbar())
  p<- ggplot(data=data, map )+plot_type+
    geom_jitter(aes_string(color=condition))+labs(main = title)+
    stat_compare_means( comparisons = comp, method = test, size=signif.size, label=label, hide.ns = hide.ns, paired=paired)

    if(!is.null(labels)){

        data[, condition]<- mgsub(string= as.character(data[,condition]),
                pattern = as.character(unique(data[,condition])),
                                   replacement = labels)## transform label as user wants it
        cat("Transformed", as.character(unique(old)), 'into -> ->', labels, "\n")#track the transformations
      v= labels
      v1=NULL
      for(i in 1:length(v)){
         v1[i]=paste0(v[i],
               ' (n=',sum(data[,condition]==v[i], na.rm = T),')')#create the labels for the y axis
      }
      p<- p+ scale_x_discrete(labels=v1)
      return(p)
    } else {
        v= factor(unique(data[,condition]))
        v= sort(v)
        x= 1:length(v)
        v1= NULL
          for(i in x) {
            if(!is.na(v[i])){
              v1[i]=paste0(v[i], ' (n=',sum(data[,condition]==v[i], na.rm = T),')')
            }
            # if(is.na(v[i])) {
            #   cat("NA here")
            #   v1[i]= paste0(v[i], " (n=", sum(is.na(data[,condition]==v[i]), na.rm = T), ')')
            # }
          }

  p <- p+scale_x_discrete(labels=v1)
      }
 return(p)
  data[,condition]<- old
  rm(old)
}#end#
