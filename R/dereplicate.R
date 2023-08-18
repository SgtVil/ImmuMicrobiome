#' Title
#'
#' @param filt_fastq_list
#' @param err_list
#' @param cores
#' @param pool
#'
#' @return
#' @export
#'
#' @examples
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
