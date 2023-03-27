#' Title
#'
#' @param df dataframe returned by \link{seq_table}
#' @param confidence_interval interval of confidence
#'
#' @return An ellipse dataframe
#' @export
#'
#' @examples No example
ellipse_me = function(df, confidence_interval){
  require(car)

  require(sp)

  ellipse = dataEllipse(df$log10_neg_abundance,
                        df$log2_neg_ratio,
                        levels=confidence_interval,
                        draw = F,
                        segments = length(df$log10_abundance)-1)
  tmp = sapply(as.list(confidence_interval), function(x){
     !point.in.polygon(df$log10_abundance,
                       df$log2_ratio,
                       ellipse[[as.character(x)]][,1],
                       ellipse[[as.character(x)]][,2])
     })
  conf = ifelse(rowSums(tmp)==3, 0.999,
                     ifelse(rowSums(tmp)==2, 0.99,
                      ifelse(rowSums(tmp)==1, 0.95, 0)))

  tmp= list(boolean = data.frame(taxonomy= df$taxonomy, ellipse_level= conf), values= as.data.frame(ellipse))
  return(tmp)
}
