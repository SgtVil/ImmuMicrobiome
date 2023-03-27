summarise_fastq<- function(fastq, cores= 1, plot=TRUE, plot_type="boxplot"){
  require(ShortRead)
  require(parallel)
  require(ggplot2)
  theme_set(theme_classic()+
              theme(plot.title = element_text(size = cores, face="bold"),
                    axis.title = element_text(size=20)))
if(length(fastq)>2){
  print("pair end analysis")
  x_fwd <- mclapply(X= fastq$fastq_fwd, readFastq, mc.cores = cores)
  x_rev = mclapply(X= fastq$fastq_rv, readFastq, mc.cores = cores)
  width_fwd <- sapply(x_fwd, ShortRead::width, simplify = T)
  width_rev <- sapply(x_rev, ShortRead::width, simplify = T)
  x1_fwd <- do.call("cbind", width_fwd) %>%
    as.data.frame %>%
    pivot_longer(cols = everything(), names_to = "fastq", values_to = "length")
  x1_rev <- do.call("cbind", width_rev) %>%
    as.data.frame %>%
    pivot_longer(cols = everything(), names_to = "fastq", values_to = "length")
  tmp_fwd <- summary(x_fwd)[,1]
  tmp_rev <- summary(x_rev)[,1]
  df_fwd <- mclapply(fastq$fastq_fwd, qa, mc.cores = cores)
  df_fwd <- sapply(df_fwd, function(x)summary(x[["readQualityScore"]]$quality))
  df_rev <- mclapply(fastq$fastq_rv, qa, mc.cores = cores)
  df_rev <- sapply(df_rev, function(x)summary(x[["readQualityScore"]]$quality))
  s_fwd = sapply(width_fwd, summary)
  s_rev <- sapply(width_rev, summary)
 df = data_frame(depth= as.numeric(c(tmp_fwd, tmp_rev)), quality= c(df_fwd["Mean",], df_rev["Mean",]), length=c(s_fwd["Mean",], s_rev["Mean",]),
                 sample_name= c(paste(fastq$names, "fwd", sep="_"), paste(fastq$names, "rev", sep="_")))
} else{
  print("single end analysis")
  x <- mclapply(X= fastq$fastq_fwd, readFastq, mc.cores = cores)
  width <- sapply(x, ShortRead::width, simplify = T)
  x1 <- do.call("cbind", width) %>%
    as.data.frame %>%
    pivot_longer(cols = everything(), names_to = "fastq", values_to = "length")
  tmp <- summary(x)[,1]
  df <- mclapply(fastq$fastq_fwd, qa, mc.cores = cores)
  df <- sapply(df, function(x)summary(x[["readQualityScore"]]$quality))
  s <- sapply(width, summary)
  df= data.frame(depth= as.numeric(tmp), quality= df["Mean",], length=s["Mean",], sample_name= fastq$names)
}

  if(plot==TRUE & plot_type=="density"){
     p <- ggplot(x1, aes(length,fill=fastq))+
      geom_density( alpha=0.2, linetype="dashed")+
      labs(y="Read count",
           title = paste(length(fastq), "files used to compute"),
           x= "Read length")
    ggsave(paste0(unique(basename(dirname(fastq))), "_length.png"), plot = p, )
  }
  if(plot==TRUE & plot_type=="boxplot"){
    p <- ggplot(x1, aes(length,fastq, fill=fastq))+
      geom_boxplot( alpha=0.2)+
      scale_x_continuous(limits = quantile(x1$length, c(0.1, 0.9)))
      labs(y="Read count",
           title = paste(length(fastq), "files used to compute"),
           x= "Read length")
    ggsave(paste0(unique(basename(dirname(fastq))), "_length.png"), plot = p )
  }


  return(df)
}
