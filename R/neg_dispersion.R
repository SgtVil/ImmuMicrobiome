#' Check the dispersion of your samples and the impact of the zero_treatment.
#'
#' Check the dispersion of your samples and the impact of the zero management you want to perform.
#' There are three possibilities [random_generation()], [no_zero()] and [minimum_count()].
#' It will produce a list of plots were you will find the \eqn{neg1/neg2} ratio and the \eqn{pos/neg1} ratio.
#' This is particularly important to check once you're working on a new set of data.
#'
#' Take into account that random_generation will be usefull for samples with a few taxa. If you samples are from healthy adults you might want to use
#'
#' @param df A element of the [seq_table()] returned object. You can run this function using a loop or lapply to make the plots for every samples.
#' @param positive_sorted_sample The pattern of the positive fraction. Default = "pos"
#' @param negative_sorted_sample The pattern of the first negative fraction. Default = "neg1"
#' @param second_negative_sample The pattern of the second negative fraction. Default = "neg2"
#' @param type Either a facet of the three possibilities as a single panel: meaning the [random_generation()],
#' the [minimum()] or [no_zero()] functions as a single panel
#' or superposed : the \eqn{pos/neg1} ratio for each way of dealing with zero values on one panel and all the \eqn{neg1/neg2} ratio on an other panel.
#'
#' @return Return a ggplot2 graph.
#' \item{facet_all_three}{For the facet_all_three : the black dots will represent the [no_zero()], the green dots the [random_generation()]
#' and the red dots the [minimum_count()].}
#' \item{superposed}{For the superposed : the grey dots represent the \eqn{neg1/neg2} and the red dots represent \eqn{pos/neg1}.}
#'
#' @export
#' @import ggpubr
#' @examples Check url to come.
neg_dispersion <- function(df, positive_sorted_sample="pos",
                           negative_sorted_sample="neg1",
                           second_negative_sample="neg2",
                           type= "facet_all_three"){

   res = list()
   pos = df[, which(grepl(colnames(df), pattern = positive_sorted_sample))]
   neg1 = df[, which(grepl(colnames(df), pattern = negative_sorted_sample))]
   neg2 = df[, which(grepl(colnames(df), pattern = second_negative_sample))]
   res$zero_ASV = sum(c( length(pos[pos==0]), length(neg1[neg1==0]), length(neg2[neg2==0]) ))

  # RANDOM GENERATION
  res$random_generation = random_generation(df, positive_sorted_sample, negative_sorted_sample, second_negative_sample)

  # MINIMUM COUNT

  res$minimum = minimum_count(df, positive_sorted_sample, negative_sorted_sample, second_negative_sample)

  # NO ZERO

  res$no_zero = no_zero(df, positive_sorted_sample, negative_sorted_sample, second_negative_sample)

  name = str_split(colnames(df)[3], pattern = "_", simplify = T)[,1]
 if(type == "facet all three"){
   p1 <- ggplot(res$random)+
     geom_jitter(aes(x=log10_neg_abundance, log2_neg_ratio), color="grey")+
     stat_ellipse(type = "norm", level = 0.95, aes(x=log10_neg_abundance, log2_neg_ratio))+
     geom_jitter(aes(x=log10_abundance, log2_ratio), color="red")+
     scale_y_continuous(limits=c(-5,5))+
     labs(title= "random generation")
   p2 <- ggplot(res$minimum )+
     geom_jitter(aes(x=log10_neg_abundance, log2_neg_ratio), color="grey")+
     stat_ellipse(type = "norm", level = 0.95, aes(x=log10_neg_abundance, log2_neg_ratio))+
     geom_jitter(aes(x=log10_abundance, log2_ratio), color="red")+
     scale_y_continuous(limits=c(-5,5))+
     labs(title= "minimum")
   p3 <- ggplot(res$no_zero )+
     geom_jitter(aes(x=log10_neg_abundance, log2_neg_ratio), color="grey")+
     stat_ellipse(type = "norm", level = 0.95, aes(x=log10_neg_abundance, log2_neg_ratio))+
     geom_jitter(aes(x=log10_abundance, log2_ratio), color="red")+
     scale_y_continuous(limits=c(-5,5))+
     labs(title= "no zero")
   facet = ggarrange(p1, p2, p3, nrow = 1 )
   p= annotate_figure(p = facet, top = text_grob(paste(name, res$zero_ASV, "ASV are equal to zero")))

 }

  if(type == "superposed"){
    p1 <- ggplot()+
      geom_jitter(data= res$random_generation , aes(x=log10_abundance, log2_ratio),
                  color="green")+
      geom_jitter(data= res$minimum ,aes(x=log10_abundance, log2_ratio), color="red")+
      geom_jitter(data = res$no_zero,aes(x=log10_abundance, log2_ratio), color="black")+
      scale_y_continuous(limits=c(-5,5))+
      labs(title= "pos on neg")
    p2 <- ggplot()+
      geom_jitter(data= res$random_generation , aes(x=log10_neg_abundance, log2_neg_ratio),
                  color="green")+
      geom_jitter(data= res$minimum ,aes(x=log10_neg_abundance, log2_neg_ratio), color="red")+
      geom_jitter(data = res$no_zero,aes(x=log10_neg_abundance, log2_neg_ratio), color="black")+
      scale_y_continuous(limits=c(-5,5))+
      labs(title= "neg on neg")
    superposed = ggarrange(p1, p2, nrow = 1 )
    p= annotate_figure(p =superposed , top = text_grob(paste(name,  res$zero_ASV, "ASV are equal to zero")))

  }
return(p)
    # return(res)

}

