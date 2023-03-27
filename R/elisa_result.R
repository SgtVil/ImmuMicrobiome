elisa_results <- function(df, absorbance, slope, dilution){
  abs = as.matrix(df[, absorbance])
  tmp=  ED(object = slope, abs[,1], display= F, type= "absolute")
  # tmp %>% filter()
  tmp = cbind(df, tmp)
  tmp$Estimate= ifelse(tmp$Estimate< slope$coefficients[2], 0, tmp$Estimate)
  tmp$tube_conc = tmp$Estimate* tmp[,dilution]
  
  rownames(tmp) = rownames(df)
  return(tmp)
}
