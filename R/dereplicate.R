#' Dereplicate sequences in a dada2 manner
#'
#' @description
#' A wrapper function for \code{\link[dada2]{dada()}}.
#'
#' BE CAREFULL, THIS FUNCTIONS HAS BEEN WRITTEN FOR LARGE AMOUNT OF RAM.
#' Currently the function might not scale efficiently for laptops.
#'
#'
#' @param filt_fastq_list List of filtered fastq returned by \link{filt_list()}
#' @param err_list Error object returned by \link{error_check()}
#' @param cores Number of threads to use. Default = 1
#' @param pool Pool, pseudo-pool or no pool. See \code{\link[dada2]{dada()}}. Default = F.
#'
#' @return
#' An list of \link{dada} objects.
#' @export
#'
#' @examples
#' # See Vignette.
dereplicate <- function(filt_fastq_list, err_list, cores=1, pool=F){
  if(!is.list(filt_fastq_list)){
   dd= dada(filt_fastq_list, err = err_list, pool = pool, multithread = cores)
   saveRDS(object = dd, file = "dd.rds")
   return(dd)
  }
  if(is.list(filt_fastq_list)){
   ddFs= dada(filt_fastq_list$fwd, err = err_list$errF, pool = pool, multithread = cores)
   saveRDS(object = ddFs, file = "ddFs.rds")
   ddRs= dada(filt_fastq_list$rev, err = err_list$errR, pool = pool, multithread = cores)
   saveRDS(object = ddRs, file = "ddRs.rds")
   return(dd=list(ddFs=ddFs, ddRs=ddRs))
  }
}
