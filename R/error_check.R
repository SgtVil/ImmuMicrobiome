#' Error estimation according to dada2
#' @description
#' A wrapper function that will directly pass argument to  \code{\link[dada2]{learnErrors}}. If `save_plot=T` the plots will
#'
#'
#'
#' @param filt_fastq_list List of filtered fastq returned by  \code{\link[ImmuMicrobiome]{filt_list}}.
#' @param cores Number of threads to use. Default = 1.
#' @param save_plot Logical. Save the plot in the current directory.
#' @param name Name for the plot saving. If NULL the plot will be named only with `Sys.Date` and "errors.png".
#'
#' @return
#' Objects from the \link{learnErrors} function. If save_plot selected then plots based based on \link{plotErrors} function will be saved
#' @export
#'
#' @examples
error_check <- function(filt_fastq_list, cores=1, save_plot=F, name=NULL){

  if(length(filt_fastq_list)!=2){
    message("single-end or merged fastq error learning")
    err= learnErrors(filt_fastq_list, verbose = T, multithread = cores)
    if(save_plot==T){
    p= plotErrors(err, err_in = T, nominalQ = T)
    if(is.null(name)) {
      ggsave(filename = paste0(Sys.Date(), "errors.png"), device = "png")}
    else {
      ggsave(filename = paste0(name, "errors.png"), device = "png")
    }
    }

    return(err)
  }
  if(length(filt_fastq_list)==2){
    message("pair-end error learning")
    errF = learnErrors(filt_fastq_list$fwd, verbose = T, multithread = cores)
    errR = learnErrors(filt_fastq_list$rev, verbose = T, multithread = cores)
    if(save_plot==T){
      p1= plotErrors(errF, err_in = T, nominalQ = T)
      p2= plotErrors(errR, err_in = T, nominalQ = T)
      p = ggpubr::ggarrange(p1, p2, nrow=1)
      if(is.null(name)) {
        ggsave(filename = paste0(Sys.Date(), "errors.png"), device = "png")}
      else {
        ggsave(filename = paste0(name, "errors.png"), device = "png")
      }
    }

    return(list(errF=errF, errRs=errR))
  }
}
