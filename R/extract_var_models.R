extract_var_models = function(obj, ncomp=2){

  if(class(tmp$model$family_id)[1]== "plsda"){
    df= obj$model
    df = sapply(df, "[", "prop_expl_var", simplify = T)
    df = sapply(df, "[", "X", simplify=T)
    df = do.call("rbind", df)
    rownames(df)= names(obj$model)
    df =  df %>%
      as.data.frame

    lim= max(apply(df, 2, max))
  }


  plot_var = function( yvar) {
      ggplot(data=df)+
      geom_col(aes_(x=fct_reorder(rownames(df), desc(df[,yvar])), y= as.name(yvar)))+
      scale_y_continuous(labels = scales::percent, limits= c(0, lim))+
      theme(axis.text.x = element_text(angle=90))
  }

  p= lapply(colnames(df)[1:ncomp], plot_var)
  p=do.call("ggarrange", c(p, ncol=ncomp))
  return(list(df, p))
}
