estimate_slope <- function(df, concentration_value, absorbance, plot=TRUE ){
  df = as.data.frame(df)
  if(!is.numeric(df[,concentration_value])){
    df[, concentration_value]= as.numeric(df[, concentration_value])
    abs = df[,absorbance] %>% as.matrix
    conc= df[, concentration_value] %>% as.matrix
    model = drm(abs~conc, fct=LL.4(names= c("slope","lower","upper", "EC50")), data = df)
  }
  abs = df[,absorbance] %>% as.matrix
  conc= df[, concentration_value] %>% as.matrix
  model = drm(abs~conc, fct=LL.4(names= c("slope","lower","upper", "EC50")), data = df)

  if(plot==TRUE){
   plot(model)
  }
  return(model)
}

# estimate_slopes <- function(df, concentration_value, absorbances){
#   if(!is.numeric(df[,concentration_value])){
#     df[, concentration_value]= as.numeric(df$concentration_value)
#   }
#   map(absorbances, drm())
# }

elisa_results <- function(df, absorbance, slope, dilution){
  abs = as.matrix(df[, absorbance])
  tmp=  ED(object = slope, abs[,1], display= F, type= "absolute")
  # tmp %>% filter()
  tmp = cbind(df, tmp)
  tmp$Estimate= ifelse(tmp$Estimate< slope$coefficients[2], 0, tmp$Estimate)
  tmp$tube_conc = tmp$Estimate* tmp[,dilution]/1000

  rownames(tmp) = rownames(df)
  return(tmp)
}


