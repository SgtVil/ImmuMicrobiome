#' Use multidimensional reduction on matrices and data.frame.
#'
#' @description
#' `plot_reduction()` stands in place of \link{beta_dispersion} for non phyloseq objects.
#' The purpose is to allow users to plot main reduction algorithm on data.frames or matrices.
#'
#' @param mat A matrix with samples in rows
#' @param clinical_data An integer specifying the columns that are discrete values.
#' @param axis_x Number of the component for X axis
#' @param axis_y Number of the component for Y axis
#' @param nf Number of component to calculate
#' @param method Method of reduction to use. Currently available : "PCA", "CA", "DCA", "tsne"
#' @param type Chose the type of plot you want, choices are : "boxplot", "pure" or "arrows". Default="boxplot".
#' @param dist The distance method to be used for PCoA, NMDS and PERMANOVA analysis, distance are directly given to \code{\link[vegan]{vegdist}}. Default = "euclidean".
#' @param stat Either permanova \code{\link[vegan]{adonis2}} with default parameter or \code{\link[vegan]{envfit}} function.
#' @param group The group to colorise and put ellipse
#' @param color_vector Color vector
#' @param legend_title Legend title
#' @param lwd Line width
#' @param conf Confidence interval. Default = 0.9.
#' @param cex Size of the text
#' @param font Type of font. Default = 2 (bold).
#' @param pch Shape
#' @param draw Categorical, draw ellipses as lines or as polygons
#' @param ylimits @param xlimits Limits for the axis
#' @param text Should the labels be printed
#' @param ncol Number of columns for the legend. Default = 1
#' @param x.intersp Adjust the legend: character interspacing factor for horizontal (x) spacing between symbol and legend text.
#'  @param y.intersp Adjust the legend: vertical (y) distances (in lines of text shared above/below each legend entry). A vector with one element for each row of the legend can be used.
#' @param where Position for the permanova value. Default = "topleft"
#' @param inset Adjust the legend: inset distance(s) from the margins as a fraction of the plot region when legend is placed by keyword.
#' @param pca Logical. For tsne only. Does the \link{tsne} need to be run on a PCA first ?
#' @param stat.cex Size of the stats text.
#' @param ...

#'
#'
#' @return A plot rendered using the `base` package and vegan functions :
#' - \code{\link[vegan]{ordispider}}
#' - \code{\link[vegan]{ordiellipse}}
#' @seealso  [plot_constrained_reduction()] [beta_dispersion()][plot_constrained_reduction()]
#' @import vegan
#' @export
#'
#' @examples
#' data(metabolomic)
#' # Reduction with boxplot on the sides
#'
#' plot_reduction(mat= metabolomic, clinical_data = 1:4, axis_x=1, axis_y=2, nf= 5,
#'  method= "PCA", type= "boxplot", group="birth_type", stat="permanova")
#'
#' # Reduction without boxplots
#'  plot_reduction(mat= metabolomic, clinical_data = 1:4, axis_x=1, axis_y=2, nf= 5,
#'  method= "PCA", type= "pure", group="birth_type", stat="permanova")
#'
#'# Reduction with arrows for the vectors loading
#'plot_reduction(mat= metabolomic, clinical_data = 1:4, axis_x=1, axis_y=2, nf= 5,
#'  method= "PCA", type= "arrows", group="birth_type", stat="permanova")
#'
#'# Test permanova or envift
#' plot_reduction(mat= metabolomic, clinical_data = 1:4, axis_x=1, axis_y=2, nf= 5,
#'  method= "PCA", type= "arrows", group="birth_type", stat="envfit")
#'
#'
#'
plot_reduction = function(mat,  clinical_data, axis_x=1, axis_y=2, nf= 5, method= "PCA", type= "boxplot",
                          group=NULL,  dist = "euclidean", stat= "none",
                          color_vector= c("cyan4","brown","deepskyblue", "black","red"),
                          legend_title= NULL, lwd=1, conf=0.9, cex=2,
                          font=2, pch=20, draw= "lines",
                          ylimits="auto", xlimits= "auto", text=F, ncol=1,
                          x.intersp = 1, y.intersp=0.5,
                          where="topleft", inset=0.2,  pca= T, scale =T,  stat.cex= 2, ...){
  old.par = par()
  on.exit(layout(matrix(c(1,1))))

  if(class(dist)=="dist"){
    d= dist
  } else {
    d= vegan::vegdist(mat[, -clinical_data], method = dist)
  }

  if(is.null(group)){
    stop("Need factor to segregate result.")
  }

  # group= as.factor(group)
  fac= mat[,group]
  fac= as.factor(fac)

if( method== "PCoA"){
  p= vegan::wcmdscale(d, eig = T)
  p_li= p$points
}

if(method=="NMDS"){
  p= metaMDS(mat[, -clinical_data], distance= dist)
  p_li= p$points
  ca1 = scores(p, display = "species", scaling="species")[,1]
  ca2 = scores(p, display = "species", scaling="species")[,2]
 }

  if(method== "PCA"){
    # p= dudi.pca(as(mat, 'matrix'), scannf = F, nf=5, scale = T, center = T)
    p = vegan::rda(as(mat[, -clinical_data], "matrix"), scale= scale)
    p_li = p$CA$u
    # p_li= p$li
    ca1 = scores(p, display = "species", scaling="species")[,1]
    ca2 = scores(p, display = "species", scaling="species")[,2]
  }

  if(method=="CA"){
    p = cca( as(mat[, -clinical_data], 'matrix'))
    p_li= p$CA$u
    ca1 = scores(p, display = "species", scaling="species")[,1]
    ca2 =scores(p, display = "species", scaling="species")[,2]
  }

  if(method=="DCA"){
    p = decorana( as(mat[, -clinical_data], 'matrix'))
    p_li= scores(p)
    ca1 = scores(p, display = "species", scaling="species")[,1]
    ca2 = scores(p, display = "species", scaling="species")[,2]
  }
  if(method=="tsne"){
    p = Rtsne::Rtsne( as(mat[, -clinical_data], 'matrix'), pca=pca)
    p_li= p$Y
  }


  if(xlimits=="auto"){
    xlimits = c(p_li[,axis_x] %>% max  /0.5, p_li[,1] %>% min  /0.8)
  }

  if(ylimits=="auto"){
    ylimits = c(p_li[,axis_y] %>% max  /0.5, p_li[,2] %>% min  /0.8)
  } else {
    ylimits= ylim
    xlimits= ylim
  }

  if(stat=="permanova"){
    res =adonis2(as.formula(as.formula(paste0("mat[, -clinical_data] ~" , group))),
                 data = as(mat, "data.frame"),
                 permutations = 999, na.action = na.exclude,
                 method = dist)
    p.val= paste("PERMANOVA\np=",res$`Pr(>F)`[1])
  }

  if(stat== "envfit") {
    res= vegan::envfit(formula=as.formula(paste0("p ~", group)), data=mat)
    p.val = paste("Goodness of fit\np=", res$factors$pvals)
  }
  col1= color_vector
  col1= col1[unique(fac)]
  col2=  color_vector
  col2= col2[fac]
  # #
  if(type=="boxplot"){
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
    disp= ordispider(p_li, groups = fac, col = adjustcolor(color_vector, alpha=0.3), lwd=lwd,  ylim=ylimits, xlim=xlimits)
    ordiellipse(p_li, groups= fac, conf= conf, col = color_vector, lwd = lwd, draw= draw,  ylim=ylimits, xlim=xlimits)
    points(p_li[,axis_x], p_li[,axis_y], bg= col2, xlab="", ylab="", las=2, pch=21, cex=cex, ylim=ylimits, xlim=xlimits, ...)

    if(text){
      text(x= unique(disp[,1:2]),  labels=unique(group), col="black", cex=cex, font=font)
    }

    if(stat=="permanova" | stat=="envfit"){

      legend(where, inset = 0.2, legend=p.val, bty="n", cex= stat.cex)
    }
    #2
    boxplot(p_li[,axis_x]~fac, data=p_li, horizontal=T, axes=F,  xlab=NULL, ylab=NULL,
            col= color_vector , xaxt="n", lwd=lwd/2,  ylim=xlimits)
    stripchart(p_li[,axis_x]~fac, data=p_li, method = "jitter", vertical=F, add=T, pch=21,
               col= "black", bg="gray", lwd=lwd/2,  ylim=xlimits)
    #3
    boxplot(p_li[,axis_y]~fac, data=p_li, axes=F,  ylab="", xlab="", col= color_vector, lwd=lwd/2,  ylim=ylimits)
    stripchart(p_li[,axis_y]~fac, data=p_li, method = "jitter", vertical=T, add=T, pch=21,
               col= "black", bg="gray", lwd=lwd/2,  ylim=ylimits)

    # #4
    plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab="", xlab="", axes=F)
    legend("center", legend=unique(fac), col = col1, title= legend_title, pch= 20,
           cex= cex, bty="n", ncol= ncol, y.intersp = y.intersp)

    # return(head(p$eig/sum(p$eig)*100,5))

  }
  if(type=="pure"){
    layout(matrix(c(1,1)))
    par(mar=c(0.5,0.5,1,0.5))
    plot(p_li[,axis_x], p_li[,axis_y], bg= col2, axes=F, xaxt="n", yaxt="n", xlab="", ylab="", las=2, pch=21, cex=cex, ylim=ylimits, xlim=xlimits, ...)

    disp= ordispider(p_li, groups = fac, col =  adjustcolor(color_vector, alpha=0.3), lwd=lwd,  ylim=ylimits, xlim=xlimits,)
    ordiellipse(p_li, groups= fac, conf= conf, col = adjustcolor(color_vector, alpha=0.3), lwd = lwd, draw= draw,  ylim=ylimits, xlim=xlimits)
    points(p_li[,axis_x], p_li[,axis_y], bg= col2,  xlab="", ylab="", las=2, pch=21, cex=cex, ylim=ylimits, xlim=xlimits, ...)

    text(x= unique(disp[,1:2]),  labels=unique(fac), col="black", cex=cex, font=font)

    # return( head(p$eig/sum(p$eig)*100, 5))
    if(stat=="permanova" | stat=="envfit"){

      legend(where, inset = 0.2, legend=p.val, bty="n", cex= stat.cex)
    }
  }

  if(type== "arrows"){
    if(method=="PCoA") stop("Can't plot loadings for a PCoA, use NMDS for that purpose")

    l.pos <-ca2# Create a vector of y axis coordinates
    lo <- which(ca2< 0) # Get the variables on the bottom half of the plot
    hi <- which(ca2> 0) # Get variables on the top half
    # Replace values in the vector
    l.pos <- replace(l.pos, lo, "1")
    l.pos <- replace(l.pos, hi, "3")

    plot(p, type="n",  axes=F, xlab="", ylab="",  bty="n")
    abline(h = 0, v = 0, col = "white", lwd = 3)
    # plot(p_li[ ,axis_x], p_li[ ,axis_y], xlim= c(min(ca1), max(ca1)), ylim = c(min(ca2), max(ca2)), bg= col2, axes=F, xlab="", ylab="", las=2, pch=21, cex=cex)
    disp= ordiellipse(p, groups= fac, conf= conf,  col = adjustcolor(color_vector, alpha=0.3),
                      lwd = lwd, draw= draw,  xlim= c(min(ca1), max(ca1)), ylim = c(min(ca2), max(ca2)))
    disp= ordispider(p, groups = fac, col =  adjustcolor(color_vector, alpha=0.3), lwd=lwd,  xlim= c(min(ca1), max(ca1)), ylim = c(min(ca2), max(ca2)))

    points(p, display = "sites",  bg= col2,  xlab="", ylab="", las=2, pch=21, cex=cex)
    # points(p_li[,axis_x], p_li[,axis_y], bg= col2,  xlab="", ylab="", las=2, pch=21, cex=cex, xlim= c(min(ca1), max(ca1)), ylim = c(min(ca2), max(ca2)), ...)
    arrows(x0=0, x1=ca1, y0= 0, y1= ca2, lwd=lwd/1.5)
    text(ca1, ca2, labels=names(ca1), col="black", pos=l.pos, lwd=lwd/1.5, cex=cex/1.5, font=2)

    if(text){
      text(x= unique(disp[,1:2]),  labels=unique(fac), col="black", cex=cex, font=font)
    }
    if(stat=="permanova" | stat=="envfit"){

      legend(where, inset = 0.2, legend=p.val, bty="n", cex= stat.cex)
    }
  }


}
