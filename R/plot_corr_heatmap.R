#' Plot a correlation heatmap
#'
#' @description
#' A heatmap of correlation coefficient based on matrices.
#'
#'
#' @param X Numeric matrix to be correlated to the Y matrix
#' @param Y Numeric matrix to be correlated to the X matrix
#' @param method Statistical method used to
#' @param cutoff Minimum correlation value to plot
#' @param cluster
#' @param mid
#' @param low
#' @param high
#' @param ratio
#' @param return_df
#'
#' @return
#' @export
#'
#' @examples
plot_corr_heatmap = function(X, Y, method= "spearman", cutoff=0,  cluster= F, mid= "white", low="lightcyan4", high="lightsalmon4", ratio=0.3, return_df= F){

  cr= Hmisc::rcorr(x = as.matrix(X), y = as.matrix(Y), type = method)
  r= reshape2::melt(cr$r[1:ncol(X),ncol(X)+1:ncol(Y)]) %>%
    mutate(p= melt(cr$P[1:ncol(X),ncol(X)+1:ncol(Y)], na.rm=F)$value)%>%
    rstatix::add_significance("p", symbols = c("****", "***", "**", "*", ""))%>%
    rstatix::adjust_pvalue("p", method = "fdr")%>%
    rstatix::add_significance("p.adj", symbols = c("****", "***", "**", "*", ""))

  if(cutoff!=0){
    var= r %>%
      filter(abs(value) > cutoff, p.adj < 0.05)
    r = r %>%
      filter(Var2 %in%   var$Var2 & Var1 %in% var$Var1)
  }


  r %>%
    group_by(Var2)%>%
    filter(p.adj<0.05)


    p= r%>%
    ggplot(aes(Var1, Var2, fill= value,
               height=sqrt(1-p),
               width=sqrt(1-p)))+
    geom_tile()+
    geom_text(aes(label=p.adj.signif), color = "black",
              fontface = "bold", size = 3)+
    # coord_fixed(ratio = ratio) +
    scale_fill_gradient2(midpoint = 0, mid = mid,
                         limits = c(-round(max(abs(r$value)),1),
                                    round(max(abs(r$value)),1)),
                         low =low, high = high, name="Correlation") +
    theme(plot.title = element_text(size=20, face="bold", hjust=0.5),
          axis.text.x = element_text(angle=30 , hjust=1, size=10, face="bold"),
          axis.text.y = element_text( size=10, face="bold"),
          axis.title = element_blank())

    if(cluster == T){
      tmp= r %>%
        select(Var1, Var2, value)%>%
        pivot_wider(values_from = value, names_from = Var2)%>%
        set_rownames(.$Var1)
      h= hclust(d = dist(t(tmp[,-1])), method = "ward.D2")
      h2= hclust(d = dist(tmp[,-1]), method = "ward.D2")

      p= p+
        scale_y_discrete(limits= colnames(tmp[,-1])[h$order] )+
        scale_x_discrete(limits= tmp$Var1[h2$order] )
    }


  if(return_df == T){

    return(list(r= cr$r, p= cr$P, long= r, clust= list(h1 = h, h2= h2), plot= p))
  } else return(p)
}

