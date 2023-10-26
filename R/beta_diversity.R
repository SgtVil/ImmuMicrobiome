#' Plot beta-diversity with numerous aspect.
#'
#'@param physeq (Required). A phyloseq object to analyse your data.
#'
#'@param distance A character string. The default is "wunifrac" from the phyloseq package.
#'@param method Computation method to display beta-diversity. Default : \code{"PCoA"}
#'@param group A character string from your metadata. Select one variable you want to segregate your data by.
#'@param color_vector A color vector by which you want plot. Make sure the colovector is sufficient for your parameters.
#'Default : c("cyan4","brown","deepskyblue", "black","red")
#'@param factor_to_plot (Optionnal). A character string you want as title.
#'@param second_factor_to_plot (Optionnal). A second factor.
#'@param lwd Line width.
#'@param clabel Label size.
#'@param cpoint Point size
#'@param addaxes Plot axis or not. Logical.
#'@param axesell Plot axesell. See [ade4::s.class].
#'@param csub Title size.
#'
#'@import vegan
#'@import phyloseq
#'@import ade4

#'
#'@export

beta_diversity= function(physeq,  dist= "wunifrac", nf= 5, method= "PCoA",
                         group=NULL, color_vector= c("cyan4","brown","deepskyblue", "black","red"),
                         factor_to_plot= NULL, second_factor_to_plot= NULL, lwd=1, cpoint=2,
                         clabel=1, addaxes= F, axesell=T, csub= 2, font=2, pch=20, species=F,
                         permanova= F, where="topleft", ...)
{
  old <- par()

  if(is.null(access(physeq, "tax_table"))){
    stop("Need a phyloseq object")}

  if(class(dist)=="dist"){
    d=dist
    } else {
      if(is.null(access(physeq, "phy_tree")) & dist=="Unifrac"){
        print("No phylogenetic tree in this phyloseq object, bray-curtis distance selected.")
        d= phyloseq::distance(physeq, method= "bray" )

      } else {
        d= phyloseq::distance(physeq, method= dist )
        # print(dist)
      }
    }
  res= vegan::adonis2(as.formula(as.formula(paste0("d ~" , group))),
                      data = as(sample_data(physeq), "data.frame"),
                      permutations = 999 )
  if(is.null(group)){
    stop("Need factor to segregate result.")
  }
  fac= sample_data(physeq)[,group]
  fac= as.factor(fac[[1]])
  par(lwd=lwd, font=font)
  # if(length(colnames(tax_table(physeq))<7))
  # {
  #   cat("pasted ASV vector. \n")
  #   taxa_names(physeq)<- paste(tax_table(physeq)[,"Genus"],
  #                              as.vector(paste("ASV_", c(1:length(tax_table(physeq)[, ncol(tax_table(physeq))])))))
  # }

  # if(length(method)==1){

  if( method== "PCoA"){
    # par(mfrow=c(1,1))
    p= dudi.pco(d, scannf = F, nf=nf)
    s.class(p$li, fac= fac, sub= paste(factor_to_plot, second_factor_to_plot),
            col = color_vector, grid = F, axesell = axesell, addaxes = addaxes, cpoint=cpoint, clabel=clabel, csub=csub, pch=pch, ...)
    # title( xlab = paste(p$eig[1]/sum(p$eig)*100, "%"))
    # cat("Inertia. \n")
    # return( head(p$eig/sum(p$eig)*100,5))

  }
  if(method=="BCA"){

    cat("Do you want PCA or PCoA based BCA ? Press [PCA] or [PCoA]. \n")
    Answer= readline(prompt="")
    if(Answer== "PCoA"){
      p= dudi.pco(d, scannf=F, nf=nf)
      b= bca(p, scannf=F, fac= fac)
       # par(mfrow=c(1,2))

      s.class(b$ls, fac= fac, sub=  paste(factor_to_plot, second_factor_to_plot),
              col = color_vector, grid = F, axesell = axesell, addaxes = addaxes, cpoint=cpoint,
              clabel=clabel, csub= csub, pch=pch, ...)

      # scatter(b, clabel=0.2, posieig="none")
      # cat("Inertia. \n")
      return( head(p$eig/sum(p$eig)*100,5))
      par(mfrow=c(1,1))
    }
    if(Answer=="PCA"){
      otu= as(otu_table(reverseASV(physeq)), 'matrix')
      colnames(otu)= paste0(tax_table(physeq)[,"Genus"], 1:length(tax_table(physeq)[,"Genus"]))
      p= dudi.pca(otu, scannf = F, nf=nf, scale = T, center = T)
      b= bca(p, scannf=F, fac= fac)
       # par(mfrow=c(1,2))

      s.class(b$ls, fac= fac, sub=  paste(factor_to_plot, second_factor_to_plot),
              col = color_vector, grid = F, axesell = axesell, addaxes = addaxes, cpoint=cpoint, clabel=clabel, csub= csub,
              add.plot = F, pch=pch, ...)
      # scatter(b, clab.row = 0.5, clab.col = 0.5, posieig = "none")
      # cat("Inertia. \n")
      # return( head(p$eig/sum(p$eig)*100,5))
      # reset.par()
      if(species==TRUE){
        # orditorp(b, display="species")
        s.arrow(p$co)
      }
    }

  }
  if(method=="NMDS"){
    # par(mfrow=c(1,1))
    otu= as(otu_table(reverseASV(physeq)), 'matrix')
    m= metaMDS(otu )
    s.class(m$points, fac=fac, sub=  paste(factor_to_plot, second_factor_to_plot),
            col = color_vector, grid = F, axesell = axesell, addaxes = addaxes,
            cpoint=cpoint, clabel=clabel, csub= csub, pch=pch, ...)
    if(species==TRUE){
      orditorp(m, display="species")
    }


  }
  if(method== "PCA"){
    par(mfrow=c(1,1))
    otu= as(otu_table(reverseASV(physeq)), 'matrix')
    p= dudi.pca(otu, scannf = F, nf=nf, scale = T, center = T)
    s.class(p$li, fac= fac, sub=  paste(factor_to_plot, second_factor_to_plot),
            col = color_vector, grid = F, axesell = axesell, addaxes = addaxes, cpoint=cpoint, clabel=clabel, csub= csub, pch=pch, ...)
    # cat("Inertia. \n")
    return( head(p$eig/sum(p$eig)*100,5))
  }

  if(permanova == T){
    # legend("topright", inset = 0.2, legend=paste("PERMANOVA:", "p-value", res$`Pr(>F)`[1]), bty="n")

    legend(where, inset = 0.2, legend=paste("PERMANOVA\n", "p=", res$`Pr(>F)`[1]), bty="n", cex = cex)
  }
  #}
   on.exit(old, add=T)
}

############################################################################

##########################################################################
