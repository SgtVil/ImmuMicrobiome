#' Plot reduction
#'
#' @param mat A matrix with samples in rows
#' @param axis_x Number of the component for X axis
#' @param axis_y Number of the component for Y axis
#' @param nf Number of component to calculate
#' @param method Method of reduction to use. Currently available : "PCA", "CA", "DCA", "tsne"
#' @param group The group to colorise and put ellipse
#' @param color_vector Color vector
#' @param permanova Logical
#' @param legend_title Legend title
#' @param title Title
#' @param lwd Line width
#' @param conf Confidence interval. Default = 0.9.
#' @param cex Size of the text
#' @param font Type of font. Default = 2 (bold).
#' @param pch Shape
#' @param draw Categorical, draw ellipses as lines or as polygons
#' @param ylimits @param xlimits Limits for the axis
#' @param text
#' @param ncol Number of columns for the legend. Default = 1
#' @param x.intersp Adjust the legend: character interspacing factor for horizontal (x)
#' spacing between symbol and legend text.
#'  @param y.intersp Adjust the legend: vertical (y) distances (in lines of text shared above/below each legend entry). A vector with one element for each row of the legend can be used.
#' @param where Position for the permanova value. Default = "topleft"
#' @param inset Adjust the legend: inset distance(s) from the margins as a fraction of the plot region when legend is placed by keyword.
#' @param boxplot Logical. Plot the marginal boxplots or not. Default= T
#' @param pca Logical. For tsne only. Does the tsne need to be run on a PCA first ?
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_reduction = function(mat, axis_x=1, axis_y=2,  nf= 5, method= "PCA",
                          group=NULL, permanova=F,
                          color_vector= c("cyan4","brown","deepskyblue", "black","red"),
                          legend_title= NULL, title=NULL, lwd=1, conf=0.9, cex=2,
                          font=2, pch=20, draw= "lines",
                          ylimits="auto", xlimits= "auto", text=F, ncol=1, species=F,
                          x.intersp = 1, y.intersp=0.5,
                          where="topleft", inset=0.2, boxplot= T, pca=T, ...){
  old.par = par()
  on.exit(layout(matrix(c(1,1))))


#
#   res= vegan::adonis2(as.formula(as.formula(paste0("d ~" , group))),
#                       data =  as(mat, 'matrix'),
#                       permutations = 999, na.action = na.exclude )

  if(is.null(group)){
    stop("Need factor to segregate result.")
  }

  group= as.factor(group)


  # if( method== "PCoA"){
  #   p= dudi.pco(d, scannf = F, nf=nf)
  #   p_li= p$li
  #
  # }
  #
  # if(method=="NMDS"){
  #   # par(mfrow=c(1,1))
  #   otu= as(otu_table(reverseASV(physeq)), 'matrix')
  #   p= metaMDS(otu)
  #   p_li= p$points
  #   if(species==TRUE){
  #     orditorp(p, display="species")
  #   }
  # }

  if(method== "PCA"){
    p= dudi.pca(as(mat, 'matrix'), scannf = F, nf=nf, scale = T, center = T)
    p_li= p$li
  }
  if(method=="CA"){
    p = cca( as(mat, 'matrix'))
    p_li= p$CA$u
  }
  if(method=="DCA"){
    p = decorana(  as(mat, 'matrix'))
    p_li= scores(p)
  }
  # if(method=="tsne"){
  #   p = Rtsne( as(mat, 'matrix'), pca=pca)
  #   p_li= p$Y
  # }

  if(xlimits=="auto"){
    xlimits = c(p_li[,axis_x] %>% max  /0.5, p_li[,1] %>% min  /0.8)
  }

  if(ylimits=="auto"){
    ylimits = c(p_li[,axis_y] %>% max  /0.5, p_li[,2] %>% min  /0.8)
  } else{
    ylimits= ylim
    xlimits= ylim
  }
  # shape= pch#set shapes
  # shape= shape[as.factor(sample_data(physeq)[,group])]
  col1= color_vector
  col1= col1[unique(group)]
  col2=  color_vector
  col2= col2[group]
  # #
  if(boxplot==T){
    # prepare the matrix
    layout(matrix(c(2,2,2,4,
                    1,1,1,3,
                    1,1,1,3,
                    1,1,1,3),
                  nrow = 4,
                  ncol = 4,
                  byrow = TRUE))

    par(mar=c(0,0,0,0))
    #1
    plot(p_li[,axis_x], p_li[,axis_y], bg= col2, axes=F, xlab="", ylab="", las=2, pch=21, cex=cex, ylim=ylimits, xlim=xlimits, ...)

    disp= ordispider(p_li, groups = group, col = adjustcolor(color_vector, alpha=0.3), lwd=lwd,  ylim=ylimits, xlim=xlimits)
    ordiellipse(p_li, groups= group, conf= conf, col = color_vector, lwd = lwd, draw= draw,  ylim=ylimits, xlim=xlimits)
    points(p_li[,axis_x], p_li[,axis_y], bg= col2, xlab="", ylab="", las=2, pch=21, cex=cex, ylim=ylimits, xlim=xlimits, ...)
    if(text){
      text(x= unique(disp[,1:2]),  labels=unique(group), col="black", cex=cex, font=font)
    }
    if(permanova == T){


      legend(where, inset = 0.2, legend=paste("PERMANOVA\n", "p=", res$`Pr(>F)`[1]), bty="n", cex=cex)
    }
    #2
    boxplot(p_li[,axis_x]~group, data=p_li, horizontal=T, axes=F,  xlab=NULL, ylab=NULL,
            col= color_vector , xaxt="n", lwd=lwd/2,  ylim=xlimits)
    stripchart(p_li[,axis_x]~group, data=p_li, method = "jitter", vertical=F, add=T, pch=21,
               col= "black", bg="gray", lwd=lwd/2,  ylim=xlimits)
    #3
    boxplot(p_li[,axis_y]~group, data=p_li, axes=F,  ylab="", xlab="", col= color_vector, lwd=lwd/2,  ylim=ylimits)
    stripchart(p_li[,axis_y]~group, data=p_li, method = "jitter", vertical=T, add=T, pch=21,
               col= "black", bg="gray", lwd=lwd/2,  ylim=ylimits)

    # #4
    plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab="", xlab="", axes=F)
    legend("center", legend=unique(group), col = col1, title= legend_title, pch= 20,
           cex= cex, bty="n", ncol= ncol, y.intersp = y.intersp)

    return(head(p$eig/sum(p$eig)*100,5))

  }
  if(boxplot==F){
    layout(matrix(c(1,1)))
    par(mar=c(0.5,0.5,1,0.5))
    #1
    plot(p_li[,axis_x], p_li[,axis_y], bg= col2, axes=F, xlab="", ylab="", las=2, pch=21, cex=cex, ylim=ylimits, xlim=xlimits, ...)

    disp= ordispider(p_li, groups = group, col = adjustcolor(color_vector, alpha=0.3), lwd=lwd,  ylim=ylimits, xlim=xlimits,)
    ordiellipse(p_li, groups= group, conf= conf, col = color_vector, lwd = lwd, draw= draw,  ylim=ylimits, xlim=xlimits)
    points(p_li[,axis_x], p_li[,axis_y], bg= col2,  xlab="", ylab="", las=2, pch=21, cex=cex, ylim=ylimits, xlim=xlimits, ...)

    text(x= unique(disp[,1:2]),  labels=unique(group), col="black", cex=cex, font=font)
    if(!is.null(title)){
      title(title)
    }
    if(permanova == T){
      # legend("topright", inset = 0.2, legend=paste("PERMANOVA:", "p-value", res$`Pr(>F)`[1]), bty="n")

      legend(where, inset = 0.2, legend=paste("PERMANOVA\n", "p=", res$`Pr(>F)`[1]), bty="n", cex=cex)
    }


    return( head(p$eig/sum(p$eig)*100, 5))

  }
}
