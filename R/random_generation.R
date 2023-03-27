random_generation = function(df, positive_sorted_sample, negative_sorted_sample, second_negative_sample){
  theme_set(theme_bw())
  pos = df[, which(grepl(colnames(df), pattern = positive_sorted_sample))]
  neg1 = df[, which(grepl(colnames(df), pattern = negative_sorted_sample))]
  neg2 = df[, which(grepl(colnames(df), pattern = second_negative_sample))]
  
  res = list()
  vec = pos!=0 & neg1!=0 & neg2!=0
  res$zero_ASV = sum(c( length(pos[pos==0]), length(neg1[neg1==0]), length(neg2[neg2==0]) ))
  
  # RANDOM GENERATION
  random= sapply(df[vec,2:4], quantile, seq(0, 1, 1/20))[1,] %>%
    mean() # create the random value
  
  pos[pos==0] <- runif(length(pos[pos==0]),
                           min= min(pos[pos!=0]),
                           max= random)
  neg1[neg1==0] <- runif(length(neg1[neg1==0]),
                             min= min(neg1[neg1!=0]),
                             max= random)
  neg2[neg2==0] <- runif(length(neg2[neg2==0]),
                             min= min(neg2[neg2!=0]),
                             max= random)
  
  log10_abundance <- log10(pos*neg1)
  log2_ratio <- log2(pos/neg1)
  log2_neg_ratio <- log2(neg1/neg2)
  log10_neg_abundance <- log10(neg1*neg2)
  
  random <- df %>%
    mutate(as.data.frame(cbind(log10_abundance, log2_ratio, log10_neg_abundance, log2_neg_ratio))) %>%
    arrange(desc(log10_abundance))
 
  return(random)
  
}
