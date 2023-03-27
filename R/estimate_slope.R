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