#' Make a constrained analysis of beta diversity.
#' @description A wrapper function to use cca, rda and capscale in the same manner as beta_dispersion.
#'
#' @param physeq A phyloseq object.
#' @param dist Either a distance mecthod given as string or a distance object.
#'   Default= "bray".
#' @param model Character vector to create a model. This is related to the functions rda, cca and capscale.
#' @param nf Number of component to keep. Default= 5.
#' @param type  Chose the type of plot you want, choices are : "boxplot", "pure" or "arrows". Default="boxplot".
#' @param stat Either permanova \code{\link[vegan]{adonis2}} with default parameter or \code{\link[vegan]{envfit}} function.
#' @param method The mathematical method to reduce the dimensions. Currently this function supports CCA, RDA and dbRDA
#' @param group String defining the groups you want to plot on the beta diversity.
#' @param color_vector Optionnal. A color vector you wish to use for plotting
#' @param legend_title Optionnal. A legend title.
#' @param conf Interval of confidence for the ellipses. Default= 0.9.
#' @param draw Either draw lines or polygon. Default = "lines".
#' @param ncol Number of columns for the legend. Default = 1.
#' @param x.intersp character interspacing factor for horizontal (x) spacing between symbol and legend text.
#' @param y.intersp vertical (y) distances (in lines of text shared above/below each legend entry). A vector with one element for each row of the legend can be used.
#' @param where Where you want the permanova result.
#' @param axis_x Which component for axis x. Default = 1
#' @param axis_y Which component for axis y. Default = 2
#' @param lwd Line width
#' @param cex Text size
#' @param font Font type
#' @param pch Shape type
#' @param ylimits Limits for y axis
#' @param xlimits Limits for x axis
#' @param text Boolean. Plot the names of the group on the graph.
#' @param species Not implemented.
#' @param inset Permanova legend adjustment see \link{legend}
#' @param ... Parameters passed to base plotting.
#'
#' @return Return a plot and the result of the constrained analysis.
#' @import phyloseq
#' @import magrittr
#' @export
#'
#' @examples
#' res= constrained_beta_dispersion(enterotype,
#'                                  model= "SeqTech+Gender+Nationality+Age+ClinicalStatus",
#'                                  group="Enterotype", method = "CCA",
#'                                  boxplot =T,
#'                                  text=T,
#'                                  color_vector = tol21rainbow)
#' #res
constrained_beta_dispersion= function(physeq, axis_x=1, axis_y=2, model, dist= "bray", nf= 5, method= "CCA",
                                      type="boxplot", stat="none",
                                      group=NULL, color_vector= c("cyan4","brown","deepskyblue", "black","red"),
                                      legend_title= NULL, lwd=1, conf=0.9, cex=2,
                                      font=2, pch=20, draw= "lines",
                                      ylimits="auto", xlimits= "auto", text=F, ncol=1, species=F,
                                      x.intersp = 1, y.intersp=0.5, where="topleft", inset=0.2,...){
  # asv= as(otu_table(reverseASV(physeq)), 'matrix')

  old.par = par()
  on.exit(layout(matrix(c(1,1))))


  mod = as.formula(paste("as(otu_table(reverseASV(physeq)), 'matrix')", "~",model))
  df= as(sample_data(physeq), "data.frame")


  fac= sample_data(physeq)[,group]
  fac= as.factor(fac[[1]])


  if(method=="CCA"){
    p= cca(mod, data=df, na.action=na.exclude)
    p_li= scores(p, display = "sites", choices=c(axis_x, axis_y) )
  }

  if(method=="RDA"){
    p= rda(mod, data=df, na.action=na.exclude)
    p_li= scores(p, display = "sites", choices=c(axis_x, axis_y) )
  }

  if(method=="dbRDA"){
    p = capscale(mod, data=df, na.action=na.exclude, dist=dist)
    p_li= scores(p, display = "sites", choices=c(axis_x, axis_y) )
  }

  if(xlimits=="auto"){
    xlimits = c(p_li[!is.na(p_li[,axis_x]),axis_x] %>% max  /0.5, p_li[!is.na(p_li[,axis_x]),axis_x] %>% min  /0.8)
  }

  if(ylimits=="auto"){
    ylimits = c(p_li[!is.na(p_li[,axis_y]),axis_y]  %>% max  /0.5, p_li[!is.na(p_li[,axis_y]),axis_y]  %>% min  /0.8)
  } else{
    ylimits= ylim
    xlimits= ylim
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
  return(p)
}
