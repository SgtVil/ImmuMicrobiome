#' A wrapper function to plot the quality of fastq.
#'
#' @param fastq_list An object given by list_fastq()
#' @param n Number of fastq you want to incorpore
#' @param save_plot Save the plot as .png plot or not. Default= FALSE
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples No example.
qc_check <- function(fastq_list, n=30, save_plot=F){

  if(length(fastq_list)<3){
    rng= sample(1:length(fastq_list[[1]]), size = n)
    p=  plotQualityProfile(fastq_list[[1]][rng], aggregate = T)
    if(save_plot==T){
    ggsave(paste0(unique(basename(dirname(fastq_list[[1]]))), "QC check.png"))
    } else {
      return(p)
    }

  }
  if(length(fastq_list)>2){
    rng= sample(1:length(fastq_list[[1]]), size = n)
    p1 = plotQualityProfile(fastq_list[[1]][rng], aggregate = T)
    p2 = plotQualityProfile(fastq_list[[2]][rng], aggregate = T)
    p= ggpubr:: ggarrange(p1, p2, nrow = 1)
    if(save_plot==T){
      ggsave(paste0(unique(basename(dirname(fastq_list[[1]]))), "QC check.png"))
    } else {
      return(p)
    }
    return(p)
  }
}
