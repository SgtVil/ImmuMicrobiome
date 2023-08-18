#' Make a constrained analysis of beta diversity.
#' @description A wrapper function to use cca, rda and capscale in the same manner as beta_dispersion.
#'
#' @param physeq A phyloseq object.
#' @param dist Either a distance mecthod given as string or a distance object.
#'   Default= "bray".
#' @param model Character vector to create a model. This is related to the functions rda, cca and capscale.
#' @param nf Number of component to keep. Default= 5.
#' @param method The mathematical method to reduce the dimensions. Currently this function supports CCA, RDA and dbRDA
#' @param group String defining the groups you want to plot on the beta diversity.
#' @param color_vector Optionnal. A color vector you wish to use for plotting
#' @param legend_title Optionnal. A legend title.
#' @param conf Interval of confidence for the ellipses. Default= 0.9.
#' @param draw Either draw lines or polygon. Default = "lines".
#' @param ncol Number of columns for the legend. Default = 1.
#' @param permanova Make the permanova analysis. Default = FALSE.
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
#' @param boxplot Draw boxplot for each axis. Default = T
#' @param inset Permanova legend adjustment see \link{legend}
#'
#' @import ade4
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
                                      group=NULL, color_vector= c("cyan4","brown","deepskyblue", "black","red"),
                                      legend_title= NULL, lwd=1, conf=0.9, cex=2,
                                      font=2, pch=20, draw= "lines",
                                      ylimits="auto", xlimits= "auto", text=F, ncol=1, species=F,
                                      x.intersp = 1, y.intersp=0.5,
                                      permanova=F, where="topleft", inset=0.2, boxplot= T,...){
  # asv= as(otu_table(reverseASV(physeq)), 'matrix')

  old.par = par()
  on.exit(layout(matrix(c(1,1))))


  mod = as.formula(paste("as(otu_table(reverseASV(physeq)), 'matrix')", "~",model))
  df= as(sample_data(physeq), "data.frame")


  fac= sample_data(physeq)[,group]
  group= as.factor(fac[[1]])


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

    disp= ordispider(p, groups = group, col = adjustcolor(color_vector, alpha=0.3), lwd=lwd,  ylim=ylimits, xlim=xlimits)
    ordiellipse(p, groups= group, conf= conf, col = color_vector, lwd = lwd, draw= draw,  ylim=ylimits, xlim=xlimits)
    points(p_li[,axis_x], p_li[,axis_y], bg= col2, xlab="", ylab="", las=2, pch=21, cex=cex, ylim=ylimits, xlim=xlimits, ...)
    text(p, dis="bp")
    if(text){
      text(x= unique(disp[,1:2]),  labels=unique(group), col="black", cex=cex, font=font)
    }
    # if(permanova == T){
    #
    #
    #   legend(where, inset = 0.2, legend=paste("PERMANOVA\n", "p=", res$`Pr(>F)`[1]), bty="n", cex=cex)
    # }
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


    return(p)
  }
  if(boxplot==F){
    layout(matrix(c(1,1)))
    par(mar=c(0,0,0,0))
    #1
    plot(p_li[,axis_x], p_li[,axis_y], bg= col2, axes=F, xlab="", ylab="", las=2, pch=21, cex=cex, ylim=ylimits, xlim=xlimits, ...)

    disp= ordispider(p, groups = group, col = adjustcolor(color_vector, alpha=0.3), lwd=lwd,  ylim=ylimits, xlim=xlimits,)
    ordiellipse(p, groups= group, conf= conf, col = color_vector, lwd = lwd, draw= draw,  ylim=ylimits, xlim=xlimits)
    points(p_li[,axis_x], p_li[,axis_y], bg= col2,  xlab="", ylab="", las=2, pch=21, cex=cex, ylim=ylimits, xlim=xlimits, ...)
    text(p, dis="bp")
    if(text){
      text(x= unique(disp[,1:2]),  labels=unique(group), col="black", cex=cex, font=font)
    }
    # if(permanova == T){
    #   # legend("topright", inset = 0.2, legend=paste("PERMANOVA:", "p-value", res$`Pr(>F)`[1]), bty="n")
    #
    #   legend(where, inset = 0.2, legend=paste("PERMANOVA\n", "p=", res$`Pr(>F)`[1]), bty="n", cex=cex)
    # }

    legend("center", legend=unique(group), col = col1, title= legend_title, pch= 20,
           cex= cex, bty="n", ncol= ncol, y.intersp = y.intersp)

    return(p)

  }

}
