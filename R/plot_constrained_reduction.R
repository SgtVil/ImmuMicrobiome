#' Make a constrained reduction on matrices and data.frame.
#'
#' @description
#' Currently the function supports only the constrained reduction around categorical data.
#' It's not possible yet to give two different Omics datasets to find correlation between them.
#'
#' @param mat A numeric matrix with samples as rows
#' @param clinical_data An integer specifying the columns that are discrete values.
#' @param axis_x Which component for axis x. Default = 1
#' @param axis_y Which component for axis y. Default = 2
#' @param model Specify here the model to be used for regression
#' @param stat Either permanova \code{\link[vegan]{adonis2}} with default parameter or \code{\link[vegan]{envfit}} function.
#' @param type  Chose the type of plot you want, choices are : "boxplot", "pure" or "arrows". Default="boxplot". Here arrows will return the vectors/factors used for the
#' constrained analysis.
#' @param nf Number of components
#' @param method The mathematical method to reduce the dimensions. Currently this function supports CCA, RDA and dbRDA
#' @param group String d
#' efining the groups you want to plot on the beta diversity.
#' @param color_vector Optionnal. A color vector you wish to use for plotting
#' @param legend_title Optionnal. A legend title.
#' @param conf Confidence interval. Default = 0.9.
#' @param draw Categorical, draw ellipses as lines or as polygons#'
#' @param lwd Line width
#' @param cex Text size
#' @param font Font type
#' @param pch Shape type
#' @param ylimits Limits for Y axis, if not specified they will be automatically set
#' @param xlimits Limits for X axis, if not specified they will be automatically set
#' @param text Plot labels or not.
#' @param ncol Number of columns for the legend. Default = 1.
#' @param x.intersp character interspacing factor for horizontal (x) spacing between symbol and legend text.
#' @param y.intersp vertical (y) distances (in lines of text shared above/below each legend entry). A vector with one element for each row of the legend can be used.
#' @param where Where you want the permanova result.
#' @param inset Permanova legend adjustment see \link{legend}
#' @param stat.cex Size of the stats text.
#' @param ...
#'
#' @import vegan
#' @export
#'
#' @return A plot rendered using the `base` package and vegan functions :
#' - \code{\link[vegan]{ordispider}}
#' - \code{\link[vegan]{ordiellipse}}
#'
#' @seealso
#' [plot_reduction()]
#' [plot_constrained_reduction()]
#' [beta_dispersion()]
#'
#' @examples
#'
#'
#'
#' data(metabolomic)
#'
#' # RDA based on metabolomic data and constrained to breastfeeding status and gender.
#'
#' res= plot_constrained_reduction( metabolomic, clinical_data = 1:5,
#'                                  stat="envfit", method = "RDA", model = "breastfeeding+sex",
#'                                  group = "birth_type", type="arrows")
#'
#' # CCA
#' res = plot_constrained_reduction(metabolomic,  clinical_data = 1:4, stat="permanova", method = "RDA", model = "breastfeeding+sex",
#'                                  group = "birth_type", type="arrows")
#'
#'
#' # dbRDA
#' res = plot_constrained_reduction(metabolomic, clinical_data = 1:4, stat="envfit", method = "CCA", model = "breastfeeding+sex",
#'                                  group = "birth_type", type="arrows")
#' anova(res, by="term")


plot_constrained_reduction= function(mat, clinical_data, axis_x=1, axis_y=2, model, nf= 5, method= "CCA", type="boxplot", scale=T,
           group=NULL, dist= "bray", stat="permanova", color_vector= c("cyan4","brown","deepskyblue", "black","red"),
           legend_title= NULL, lwd=1, conf=0.9, cex=2,
           font=2, pch=20, draw= "lines",
           ylimits="auto", xlimits= "auto", text=F, ncol=1,
           x.intersp = 1, y.intersp=0.5, where="topleft", inset=0.2, stat.cex=2, axis.prop=T, cex.arrows=1, by="margin", ...){
    # asv= as(otu_table(reverseASV(physeq)), 'matrix')

    old.par = par()
    on.exit(layout(matrix(c(1,1))))


    mod = as.formula(paste("as(mat[, -clinical_data], 'matrix') ~", model))
    # df= as(sample_data(physeq), "data.frame")

    fac= as.factor(mat[,group])


    if(method=="CCA"){
      p= cca(mod, data=mat[,clinical_data], na.action=na.exclude, scale= scale)
      p_li= scores(p, display = "sites", choices=c(axis_x, axis_y) )

      ca1 = scores(p, display = "species", scaling="species")[,1]
      ca2 = scores(p, display = "species", scaling="species")[,2]
    }

    if(method=="RDA"){
      p= rda(mod, data=mat[,clinical_data], na.action=na.exclude, scale= scale)
      p_li= scores(p, display = "sites", choices=c(axis_x, axis_y) )

      ca1 = scores(p, display = "species", scaling="species")[,1]
      ca2 = scores(p, display = "species", scaling="species")[,2]
    }

    if(method=="dbRDA"){
      p = dbrda(mod, data=mat, na.action=na.exclude, dist="bray")
      #  p = capscale(mod, data=mat, na.action=na.exclude, dist="bray")
      p_li= scores(p, display = "sites", choices=c(axis_x, axis_y) )

    ca1 = scores(p, display = "species", scaling="species")[,1]
    ca2 = scores(p, display = "species", scaling="species")[,2]
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

    if(stat=="permanova"){
      res =adonis2(mod,
                   data = as(mat, "data.frame"),
                   permutations = 999, na.action = na.exclude,
                   method = dist, by=by)
      p.val= paste("PERMANOVA\np=",res$`Pr(>F)`[1])
    }

    if(stat== "envfit") {
      res= vegan::envfit(formula=as.formula(paste0("p ~", group)), data=mat)
      p.val = paste("Goodness of fit\np=", res$factors$pvals)
    #
    #   if(is.null(res$factor)){
    #     arr= scores(res, dis="vectors")
    #     } else{
    #       arr= scores(res, dis="factors")
    #         }
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
      # arrows(x1 = arr[,1], y1=arr[,2], x0 = 0, y0 = 0)
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

      legend(where, legend=unique(fac), col = col1, title= legend_title, pch= 20,
             cex= cex, bty="n", ncol= ncol, y.intersp = y.intersp)
      # return( head(p$eig/sum(p$eig)*100, 5))
      if(stat=="permanova" | stat=="envfit"){

        legend(where, inset = 0.2, legend=p.val, bty="n", cex= stat.cex)
      }
    }

    if(type== "arrows"){
      if(method=="PCoA") stop("Can't plot loadings for a PCoA, use NMDS for that purpose")

      plot(p_li[,axis_x], p_li[,axis_y], type="n", axes=F, xlab="", ylab="",  bty="n", xaxt="n", yaxt="n")
      abline(h = 0, v = 0, col = "white", lwd = 3)
      # plot(p_li[ ,axis_x], p_li[ ,axis_y], xlim= c(min(ca1), max(ca1)), ylim = c(min(ca2), max(ca2)), bg= col2, axes=F, xlab="", ylab="", las=2, pch=21, cex=cex)
      disp= ordiellipse(p, groups= fac, conf= conf,  col = adjustcolor(color_vector, alpha=0.3),
                        lwd = lwd, draw= draw,  xlim= c(min(ca1), max(ca1)), ylim = c(min(ca2), max(ca2)))
      disp= ordispider(p, groups = fac, col =  adjustcolor(color_vector, alpha=0.3), lwd=lwd,  xlim= c(min(ca1), max(ca1)), ylim = c(min(ca2), max(ca2)))

      points(p, display = "sites",  bg= col2,  xlab="", ylab="", las=2, pch=21, cex=cex)
      # points(p_li[,axis_x], p_li[,axis_y], bg= col2,  xlab="", ylab="", las=2, pch=21, cex=cex, xlim= c(min(ca1), max(ca1)), ylim = c(min(ca2), max(ca2)), ...)
      # arrows(x0=0, x1=ca1, y0= 0, y1= ca2, lwd=lwd/1.5)
      # text(ca1, ca2, labels=names(ca1), col="black", pos=l.pos, lwd=lwd/1.5, cex=cex/1.5, font=2)
      text(p, dis="reg", choices = c(axis_x, axis_y), cex=cex.arrows, lwd= cex.arrows*2, font=font)
      legend(where, legend=unique(fac), col = col1, title= legend_title, pch= 20,
             cex= cex, bty="n", ncol= ncol, y.intersp = y.intersp)

      if(axis.prop==T){
        eig= as.matrix(eigenvals(p)/sum(eigenvals(p)))
        title(xlab=eig[1], ylab = round(eig[2],3))
      }
      if(text){
        text(x= unique(disp[,1:2]),  labels=unique(fac), col="black", cex=cex, font=font)
      }
      if(stat=="permanova" | stat=="envfit"){

        legend(where, inset = 0.2, legend=p.val, bty="n", cex= stat.cex)
      }
    }

return(list(model = p, stat= res))
  }
