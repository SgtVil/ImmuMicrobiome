#' Plot the distribution of your data.
#'
#' @description
#' A barplot showing the mean abundance and the variance for each features.
#'
#'
#' @param mat A matrix or data.frame with samples as rows
#' @param top The top N features to keep. Default = 50
#' @param selection_var Which variable is to be used for selecting the top features. Default = "mean.abundance"
#' @param col Color vector
#' @param return_df Logical. Return the data.frame used to plot. Default = T.
#'
#' @return
#' @export
#'
#' @examples
plot_distribution = function(mat, top=50, selection_var= "mean.abun", col=c("orange", "brown"), return_df= F){

  mean.feat =data.frame(mean.abun= apply(as.matrix(mat), 2, mean),
                        variance = apply(as.matrix(mat), 2, var)) %>%
    as.data.frame %>%
    rownames_to_column()%>%
    set_colnames(c("features", "mean.abun", "variance")) %>%
    mutate(mean.abun= scale(mean.abun, center=F),
           variance=scale(variance, center=F))

 p1=  mean.feat%>%
    # inner_join(kegg, by=c("features"= "Query"))%>%
    top_n(wt = mean.feat[,selection_var], 50)%>%
    ggplot()+
    geom_col(aes(mean.abun, fct_reorder(features, mean.abun)),col="grey2", fill=col[1])+
   labs(x= "Scaled mean abundance")+
   theme_bw()+
   theme(axis.text.y = element_blank(),
         axis.ticks.y= element_blank(),
         axis.title.y = element_blank(),
         axis.title.x = element_text(size=20, face="bold"),
         panel.grid = element_blank(),
         axis.text.x = element_text(size=15, face="bold"))

 p2= mean.feat %>%
   top_n(wt =  mean.feat[,selection_var], 50)%>%
   ggplot()+
   geom_col(aes(variance,fct_reorder(features, mean.abun)),col="grey2", fill=col[2])+
   scale_x_reverse()+
   labs(x= "Scaled variance")+
   theme_bw()+
   theme(axis.title.y = element_blank(),
         axis.text.y = element_text(size=15, face="bold"),
         axis.text.x = element_text(size=15, face="bold"),
         axis.title.x = element_text(size=20, face="bold"),
         panel.grid = element_blank())
p= ggarrange(p2, ggplot() + theme_void(), p1, widths = c(3, 0, 2), nrow=1 )

if(return_df){
  return(list(plot= ggarrange(p2, ggplot() + theme_void(), p1, widths = c(3, 0, 2), nrow=1 ),
              df= mean.feat))
}
return(p)
}
