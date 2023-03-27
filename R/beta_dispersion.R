#' This a wrapper function to plot beta diversity directly from a phyloseq
#' object. This is based on the base package and mostly [vegan].
#'
#' @param physeq A phyloseq object.
#' @param dist Either a distance mecthod given as string or a distance object.
#'   Default= "bray".
#' @param nf Number of component to keep. Default= 5.
#' @param method Method use to plot the beta diversity. Currently this function
#'   supports only "PCoA".
#' @param group String defining the groups you want to plot on the beta diversity.
#' @param color_vector Optionnal. A color vector you wish to use for plotting
#' @param legend_title Optionnal. A legend title.
#' @param lwd,cex,font,pch,ylimits,xlimits Parameters imported from base plotting
#' @param conf Interval of confidence for the ellipses. Default= 0.9.
#' @param draw Either draw lines or polygon. Default = "lines".
#' @param textOptionnal. Labels to use for plotting. Default will be the name of the groups you gave.
#' @param ncol Number of columns for the legend. Default = 1.
#' @param permanova Make the permanova analysis. Default = FALSE.
#' @param x.intersp character interspacing factor for horizontal (x) spacing between symbol and legend text.
#' @param y.intersp vertical (y) distances (in lines of text shared above/below each legend entry). A vector with one element for each row of the legend can be used.
#' @param where Where you want the permanova result.
#' @param inset Permanova legend adjustment see \link{legend}
#' @param ... Passed to other packages.
#' @import ade4
#' @import phyloseq
#' @import magrittr
#' @export


beta_dispersion = function(physeq, axis_x=1, axis_y=2,  dist= "bray", nf= 5, method= "PCoA",
                           group=NULL, color_vector= c("cyan4","brown","deepskyblue", "black","red"),
                           legend_title= NULL, lwd=1, conf=0.9, cex=2,
                           font=2, pch=20, draw= "lines",
                           ylimits="auto", xlimits= "auto", text=F, ncol=1, species=F,
                           x.intersp = 1, y.intersp=0.5,
                           permanova=F, where="topleft", inset=0.2, ...){
  old.par = par()
  on.exit(layout(matrix(c(1,1))))
  # prepare the matrix
  layout(matrix(c(2,2,2,4,
                  1,1,1,3,
                  1,1,1,3,
                  1,1,1,3),
                nrow = 4,
                ncol = 4,
                byrow = TRUE))


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
  res= vegan::adonis2(as.formula(as.formula(paste0("d ~" , group))),
                      data = as(sample_data(physeq), "data.frame"),
                      permutations = 999 )

  if(is.null(group)){
    stop("Need factor to segregate result.")
  }
  fac= sample_data(physeq)[,group]
  group= as.factor(fac[[1]])
  # par(lwd=lwd, font=font)

  if( method== "PCoA"){
    # par(mfrow=c(1,1))
    p= dudi.pco(d, scannf = F, nf=nf)
    p_li= p$li
    # cat("Inertia. \n")
    # return( head(p$eig/sum(p$eig)*100,5))

  }
  if(method=="BCA"){
    cat("Do you want PCA or PCoA based BCA ? Press [PCA] or [PCoA]. \n")
    Answer= readline(prompt="")
    if(Answer== "PCoA"){
      p= dudi.pco(d, scannf=F, nf=nf)
      p= bca(p, scannf=F, fac= group)
      p_li= p$ls
      # cat("Inertia. \n")
      # return( head(p$eig/sum(p$eig)*100,5))
      # par(mfrow=c(1,1))
    }

    if(Answer=="PCA"){
      otu= as(otu_table(reverseASV(physeq)), 'matrix')
      p= dudi.pca(otu, scannf = F, nf=nf, scale = T, center = T)
      p= bca(p, scannf=F, fac= group)
      p_li= p$ls
      # par(mfrow=c(1,2))
      # cat("Inertia. \n")
      # return( head(p$eig/sum(p$eig)*100,5))
      # reset.par()
      if(species==TRUE){
        orditorp(b, display="species")
      }
    }
  }

  if(method=="NMDS"){
    # par(mfrow=c(1,1))
    otu= as(otu_table(reverseASV(physeq)), 'matrix')
    p= metaMDS(otu)
    p_li= p$points
    if(species==TRUE){
      orditorp(p, display="species")
    }
  }

  if(method== "PCA"){
    otu= as(otu_table(reverseASV(physeq)), 'matrix')
    p= dudi.pca(otu, scannf = F, nf=nf, scale = T, center = T)
    p_li= p$li
  }

  if(xlimits=="auto"){
    xlimits = c(p$li[,axis_x] %>% max  /0.8, p$li[,1] %>% min  /0.8)
  }

  if(ylimits=="auto"){
    ylimits = c(p$li[,axis_y] %>% max  /0.8, p$li[,2] %>% min  /0.8)
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

  par(mar=c(0,0,0,0))
  #1
  plot(p_li[,axis_x], p_li[,axis_y], bg= col2, axes=F, xlab="", ylab="", las=2, pch=21, cex=cex, ylim=ylimits, xlim=xlimits, ...)
  disp= ordispider(p, groups = group, col = adjustcolor(color_vector, alpha=0.3), lwd=lwd/2,  ylim=ylimits, xlim=xlimits,)
  ordiellipse(p, groups= group, conf= conf, col = color_vector, lwd = lwd, draw= draw,  ylim=ylimits, xlim=xlimits)

  if(text){
    text(x= unique(disp[,1:2]),  labels=unique(group), col="black", cex=cex, font=font)
  }
  if(permanova == T){
    # legend("topright", inset = 0.2, legend=paste("PERMANOVA:", "p-value", res$`Pr(>F)`[1]), bty="n")

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

  return( head(p$eig/sum(p$eig)*100,5))


}


