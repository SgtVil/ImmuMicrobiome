#' Analyse the ratio between your IgA positive fraction and the two IgA negative
#' fractions.
#'
#' @description  The \code{slide_z} function will make a Z test based on the
#'   log2 ratio of \eqn{pos*neg1} (\code{\link{slide_z_standard}}) and \eqn{neg1/neg2} (\code{\link{slid_z_modern}}).
#' @param df (Required). An object of class dataframe or matrix generated using
#'   \code{seq_table}
#' @param positive_sorted_sample Character present in every positive sorted
#'   samples, for example "pos"
#' @param negative_sorted_sample Character present in every negative sorted
#'   samples, for example "neg1"
#' @param second_negative_sample Character present in every negative sorted
#'   samples, for example "neg2"
#' @param zero_treatment (Default = "random generation") The treatment you wish
#'   to apply to ASV equal to zero in one or two samples in order to compare it
#'   with a non zero sample. Random generation will generate random numbers
#'   between the minimum value different from zero and the 5\% quartile of the
#'   data.
#' @param deltaX (Default = 50) The number of ASV to take at a time to create
#'   the main algorithm. See Details for more information
#' @param alpha The alpha
#' @param slide_version (Soon to be defunct) \code{\link{slide_z_modern}} or
#'   \code{\link{slide_z_standard}}. The standard version (not recommended) will
#'   use the log2(pos/neg1) ratio to create the normal distribution and
#'   deviation while slide_z_modern will use the log2(neg1/neg2) ratio to
#'   perform the same step.
#' @param confidence_interval (Default= c("0.95", "0.99", "0.999")) The
#'   confidence interval for the ellipse to return.
#' @param plot (Default= FALSE) Wether you want to print the slide_z plots. If
#'   TRUE, a plot of log2(pos/neg1) and log10(pos*neg1) with confidence ellipse
#'   will be drawn using ggplot2 and plotly with the taxonomy as tool. Plotting
#'   will slow drastically the code.
#'
#' @return Return an IgaSeq object containing : The function will return a
#'   single IgASeq for each samples in this present version of the package. You
#'   need to perform a \code{for} loop or \code{lapply} to make an IgASeq for
#'   each samples. The IgASeq object is composed of the following :
#'   \item{ig_seq_all}{All the taxa found in the phyloseq object.}
#'
#'   \item{ig_up}{All the taxa that are enriched in the positive fraction}
#'
#'   \item{ig_down}{All the taxa that are enriched in the negative fraction}
#'
#'   Additionally, each slot contains : \item{score}{The IgASeq score obtained
#'   with the \code{\link{slide_z}}}
#'   \item{ellipse}{The value of confidence for
#'   a given taxa.}
#'   \item{log_ratio values}{Refers to \code{\link{log_ratio}}
#'   for more information on the values given by log_ratio.}
#' @export
slide_z = function(df, positive_sorted_sample= "pos", negative_sorted_sample="neg1",
                   second_negative_sample= "neg2", zero_treatment="random generation",
                   deltaX=50, alpha=0.05, slide_version= "slide_z_modern",
                   confidence_interval=c(0.95,0.99, 0.999), plot=FALSE){
  require(htmlwidgets)
  require(plotly)
  require(tidyverse)
  if(dim(df)[2]>6){
    warning(paste(unique(df[,"sample_id"]), "duplicated names here, samples ignored"))
  } else {


    df = log_ratio(df, positive_sorted_sample, negative_sorted_sample, second_negative_sample,
                   zero_treatment)
    # if(dim(Backup)[2]<5 ) warning(paste(colnames(df),"There is multiple occurrence of positive and negative references in your data"))
    if(slide_version== "slide_z_standard"){
      res <-  slide_z_standard(df, log2_ratio, log10_abundance, positive_sorted_sample, negative_sorted_sample, deltaX, alpha=0.05)
    }

    if(slide_version=="slide_z_modern"){
      res <- slide_z_modern(df, log2_ratio, log10_abundance, deltaX, alpha=0.05)
    }
    ellipse = ellipse_me(df, confidence_interval)


    if(plot==TRUE){
      df= res@ig_seq_all
      df = mutate(df, tooltip = paste(taxonomy, "\nslide_z: ",
                                            round(score, digits=3)))
      ell= res@ellipse_data
      ell= ell %>% mutate(tooltip = df$tooltip)
      stat = df$score>=2 | df$score <= -2
      dir.create("./IgAseq plots/png", recursive = T)
      dir.create("./IgAseq plots/plotly", recursive = T)
      p <- ggplot(df)+
        geom_polygon(data = ell, aes(ell[,1], ell[,2]), alpha=0.2, fill="burlywood", linetype=2)+
        geom_path(data = ell, aes(ell[,3], ell[,4]), color="chocolate1", linetype=2)+
        geom_path(data = ell, aes(ell[,5], ell[,6]), color="coral3", linetype=2)+
        geom_point(aes(log10_abundance, log2_ratio, size=ellipse_level, color=as.character(ellipse_level), text=tooltip, shape=stat), alpha=0.6) +
        geom_point(aes(log10_neg_abundance, log2_neg_ratio), color="black")+
        scale_color_manual(values=c("cyan4","burlywood","chocolate1", "coral3"), name= "Confidence interval")+
        labs(x= "Log10 pos * neg", y="Log2 pos/neg")+
        theme_minimal()

      ggsave(filename = paste0("./IgAseq plots/png/", unique(df[,"sample_id"]), "_",zero_treatment ,'_.png'), plot = p)
      pp = plotly::ggplotly(p, tooltip ="text")
      saveWidget(pp, file=paste0("./IgAseq plots/plotly/",unique(df[,"sample_id"]) ,colnames(df)[2], "_vs_", colnames(df)[3], "_vs_", colnames(df)[4],"_", zero_treatment,"_ratio_intensity_ellipse.html"))
    }

    res <- ig_seq(ig_seq_all = merge(res$ig_seq_all, ellipse$boolean, by="taxonomy"),
                  ig_up = merge(res$ig_up, ellipse$boolean, by="taxonomy"),
                  ig_down = merge(res$ig_down, ellipse$boolean, by="taxonomy"),
                  ellipse_data= ellipse$values)
    return(res)

  }


}


