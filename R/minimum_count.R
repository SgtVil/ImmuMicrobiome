minimum_count = function(df, positive_sorted_sample, negative_sorted_sample, second_negative_sample){
  pos = df[, which(grepl(colnames(df), pattern = positive_sorted_sample))]
  neg1 = df[, which(grepl(colnames(df), pattern = negative_sorted_sample))]
  neg2 = df[, which(grepl(colnames(df), pattern = second_negative_sample))]
  
  pos[pos==0] <- min(pos[pos!=0])
  neg1[neg1==0] <- min(neg1[neg1!=0])
  neg2[neg2==0] <- min(neg2[neg2!=0])
  
  
  log10_abundance <- log10(pos*neg1)
  log2_ratio <- log2(pos/neg1)
  log2_neg_ratio <- log2(neg1/neg2)
  log10_neg_abundance <- log10(neg1*neg2)
  df= df %>%
    mutate(as.data.frame(cbind(log10_abundance, log2_ratio, log10_neg_abundance, log2_neg_ratio))) %>%
    arrange(desc(log10_abundance))
  return(df)
}
