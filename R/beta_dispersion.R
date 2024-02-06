#'  Make an unconstrained analysis of beta diversity.
#'
#' @description This a wrapper function to plot beta diversity directly from a phyloseq
#' object. This is based on the base package and mostly [vegan].
#'
#' @param physeq A phyloseq object.
#' @param dist Either a distance mecthod given as string or a distance object.
#'   Default= "bray".
#' @param nf Number of component to keep. Default= 5.
#' @param method Method use to plot the beta diversity. Currently this function
#'   supports only "PCoA".
#' @param group String defining the groups you want to plot on the beta diversity
#' @param stat Either permanova \code{\link[vegan]{adonis2}} with default parameter or \code{\link[vegan]{envfit}} function.
#' @param type  Chose the type of plot you want, choices are : "boxplot", "pure" or "arrows". Default="boxplot".
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
#' @param ylimits Limits for Y axis, if not specified they will be automatically set
#' @param xlimits Limits for X axis, if not specified they will be automatically set
#' @param text Should the labels be printed
#' @param species NMDS only, will plot the species loadings on the graph.
#' @param inset Permanova legend adjustment see \link{legend}
#' @param pca Argument for t-SNE, does the t-SNE needs to run on a PCA first or not ? Default=FAL
#' @param stat.cex ize of the stats text.
#' @param legend.cex Size of the legend text.
#' @param ...
#'
#' @import phyloseq
#' @import magrittr
#' @import Rtsne
#' @export
#'
#' @return A plot rendered using the `base` package and vegan functions :
#' - \code{\link[vegan]{ordispider}}
#' - \code{\link[vegan]{ordiellipse}}
#' @seealso
#'  [plot_reduction()]
#'  [plot_constrained_reduction()]
#'  [constrained_beta_dispersion()]
#' @examples
#'
#'data(enterotype)
#'# PCoA based on
#' beta_dispersion(enterotype, dist = "bray", method = "PCoA", group = "SeqTech",
#' color_vector = c("#777711", "#117777", "#DD7788"),
#' legend_title = "Sequencing tech", lwd = 2,
#' font = 2, draw = "polygon", text = T, permanova = T,
#' y.intersp = 0.7)
#'

beta_dispersion = function(physeq, axis_x=1, axis_y=2,  dist= "bray", nf= 5, method= "PCoA", type= "boxplot",
                           group=NULL,  stat= "none", color_vector= c("cyan4","brown","deepskyblue", "black","red"),
                           legend_title= NULL, lwd=1, conf=0.9, cex=2,
                           font=2, pch=20, draw= "lines",
                           ylimits="auto", xlimits= "auto", text=F, ncol=1, species=F,
                           x.intersp = 1, y.intersp=0.5, where="topleft", inset=0.2, pca=T, stat.cex= 2, legend.cex=2,...){
  old.par = par()
  on.exit(layout(matrix(c(1,1))))



  if(is.null(access(physeq, "tax_table"))){
    stop("Need a phyloseq object")}

  if(class(dist)=="dist"){
    d= dist
  } else {
    if(is.null(access(physeq, "phy_tree")) & dist=="Unifrac"){
      print("No phylogenetic tree in this phyloseq object, bray-curtis distance selected.")
      d= phyloseq::distance(physeq, method= "bray")
    } else {
      d= phyloseq::distance(physeq, method= dist)
      }
  }

  if(stat=="permanova"){
    res = adonis2(as.formula(paste0(" as(otu_table(reverseASV(physeq)), 'matrix') ~" , group)),
                 data = as( sample_data(physeq), "data.frame"),
                 permutations = 999, na.action = na.exclude,
                 method = dist)
    p.val= paste("PERMANOVA\np=",res$`Pr(>F)`[1])
  }

  if(stat== "envfit") {
    res= vegan::envfit(formula=as.formula(paste0(" as(otu_table(reverseASV(physeq)), 'matrix') ~" , group)),
                       data=as(sample_data(physeq), "data.frame"))
    p.val = paste("Goodness of fit\np=", res$factors$pvals)
  }
  if(is.null(group)){
    stop("Need factor to segregate result.")
  }
  fac= sample_data(physeq)[,group]
  fac= as.factor(fac[[1]])


  if( method== "PCoA"){
    p= vegan::wcmdscale(d, eig = T)
    p_li= p$points
  }

  if(method=="NMDS"){
    # par(mfrow=c(1,1))
    otu= as(otu_table(reverseASV(physeq)), 'matrix')
    p= metaMDS(otu)
    p_li= p$points
    ca1 = scores(p, display = "species", scaling="species")[,1]
    ca2 = scores(p, display = "species", scaling="species")[,2]
    if(species==TRUE){
      orditorp(p, display="species")
    }
  }

  if(method== "PCA"){
    otu= as(otu_table(reverseASV(physeq)), 'matrix')
    p= vegan::rda(otu, scale=T)
    p_li = p$CA$u
    ca1 = scores(p, display = "species", scaling="species")[,1]
    ca2 = scores(p, display = "species", scaling="species")[,2]
  }
  if(method=="CA"){
    p = cca( as(otu_table(reverseASV(physeq)), 'matrix'))
    p_li= p$CA$u
    ca1 = scores(p, display = "species", scaling="species")[,1]
    ca2 = scores(p, display = "species", scaling="species")[,2]
  }
  if(method=="DCA"){
    p = decorana( as(otu_table(reverseASV(physeq)), 'matrix'))
    p_li= scores(p)
    ca1 = scores(p, display = "species", scaling="species")[,1]
    ca2 = scores(p, display = "species", scaling="species")[,2]
  }
  if(method=="tsne"){
    p = Rtsne(as(otu_table(reverseASV(physeq)), 'matrix'), pca=pca)
    p_li= p$Y
  }

  if(xlimits=="auto"){
    xlimits = c(p_li[,axis_x] %>% max  /0.5, p_li[,1] %>% min  /0.8)
  }

  if(ylimits=="auto"){
    ylimits = c(p_li[,axis_y] %>% max  /0.5, p_li[,2] %>% min  /0.8)
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
           cex= legend.cex, bty="n", ncol= ncol, y.intersp = y.intersp)

    # return(head(p$eig/sum(p$eig)*100,5))

  }
  if(type=="pure"){
    layout(matrix(c(1,1)))
    par(mar=c(0.5,0.5,1,0.5))
    plot(p_li[,axis_x], p_li[,axis_y], bg= col2, axes=F, xaxt="n", yaxt="n", xlab="", ylab="", las=2, pch=21, cex=cex, ylim=ylimits, xlim=xlimits, ...)

    disp= ordispider(p_li, groups = fac, col =  adjustcolor(color_vector, alpha=0.3), lwd=lwd,  ylim=ylimits, xlim=xlimits,)
    ordiellipse(p_li, groups= fac, conf= conf, col = adjustcolor(color_vector, alpha=0.3), lwd = lwd, draw= draw,  ylim=ylimits, xlim=xlimits)
    points(p_li[,axis_x], p_li[,axis_y], bg= col2,  xlab="", ylab="", las=2, pch=21, cex=cex, ylim=ylimits, xlim=xlimits, ...)

    # text(x= unique(disp[,1:2]),  labels=unique(fac), col="black", cex=cex, font=font)
    legend(where, legend = unique(fac), col = col1, title = legend_title,
           pch = 20, cex = legend.cex, bty = "n", ncol = ncol, y.intersp = y.intersp)
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
    # text(ca1, ca2, labels=names(ca1), col="black", pos=l.pos, lwd=lwd/1.5, cex=cex/1.5, font=2)
    legend(where, legend = unique(fac), col = col1, title = legend_title,
           pch = 20, cex = legend.cex, bty = "n", ncol = ncol, y.intersp = y.intersp)

    if(text){
      text(x= unique(disp[,1:2]),  labels=unique(fac), col="black", cex=cex, font=font)
    }
    if(stat=="permanova" | stat=="envfit"){

      legend(where, inset = 0.2, legend=p.val, bty="n", cex= stat.cex)
    }
  }

}


