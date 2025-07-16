#' Plot the distribution of your data.
#'
#' @description
#' A barplot showing the mean abundance and the variance for each features.
#'
#'
#' @param mat A matrix or data.frame with samples as rows
#' @param top The top N features to keep. Default = 50
#' @param selection_var Which variable is to be used for selecting the top features. Default = "mean.abun"
#' @param col Color vector
#' @param return_df Logical. Return the data.frame used to plot. Default = T.
#'
#' @return Return a plot if `return_df=F` is not selected, return both if `return_df=T`
#' @export
#'
#' @examples TBD
#'
plot_distribution = function(mat, top=50, selection_var= "mean.abun", col=c("orange", "brown"), return_df= F){

  mean.feat =data.frame(mean.abun= apply(as.matrix(mat), 2, mean),
                        variance = apply(as.matrix(mat), 2, var)) %>%
    as.data.frame %>%
    rownames_to_column()%>%
    magrittr::set_colnames(c("features", "mean.abun", "variance")) %>%
    mutate(mean.abun= scale(mean.abun, center=F),
           variance=scale(variance, center=F))

 p1=  mean.feat%>%
    # inner_join(kegg, by=c("features"= "Query"))%>%
    top_n(wt = mean.feat[,selection_var], top)%>%
    ggplot()+
    geom_col(aes(mean.abun, fct_reorder(features, mean.abun)),col="grey2", fill=col[1])+
   labs(x= "Scaled mean abundance")+
   theme_bw()+
   theme(axis.title.y = element_blank(),
         panel.grid.major.y  = element_blank())

 p2= mean.feat %>%
   top_n(wt =  mean.feat[,selection_var], top)%>%
   ggplot()+
   geom_col(aes(variance,fct_reorder(features, mean.abun)),col="grey2", fill=col[2])+
   scale_x_reverse()+
   labs(x= "Scaled variance")+
   theme_bw()+
   theme(axis.title.y = element_blank(),
         panel.grid.major.y = element_blank())

p= ggarrange(p2, ggplot() + theme_void(), p1, widths = c(3, 0, 2), nrow=1 )

if(return_df){
  return(list(plot=p,
              df= mean.feat))
}
return(p)
}
