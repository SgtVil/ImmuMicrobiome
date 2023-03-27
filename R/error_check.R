error_check <- function(filt_fastq_list, cores=1, save_plot=F, name=NULL){

  if(length(filt_fastq_list)!=2){
    message("single-end or merged fastq error learning")
    err= learnErrors(filt_fastq_list, verbose = T, multithread = cores)
    if(save_plot==T){
    p= plotErrors(err, err_in = T, nominalQ = T)
    if(is.null(name)) {
      ggsave(filename = paste0(dirname(name), "errors.png"), device = "png")}
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
        ggsave(filename = paste0(dirname(name), "errors.png"), device = "png")}
      else {
        ggsave(filename = paste0(name, "errors.png"), device = "png")
      }
    }

    return(list(errF=errF, errRs=errR))
  }
}
