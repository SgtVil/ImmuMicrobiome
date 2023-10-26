#' Title
#'
#' @param fastq_list List of fastq returned by \link{list_fastq}.
#' @param filt_list List of filtered fastq names givent by \link{filt_list}.
#' @param cutting_param Equivalent to `truncLen` in \link{filterAndTrim}.
#' @param maxEE Equivalent to `maxEE` in \link{filterAndTrim}.
#' @param trimleft
#' @param cores
#'
#' @return
#' @export
#'
#' @examples
filter_fastq = function(fastq_list, filt_list,  cutting_param, maxEE, trimleft=0, cores,... )
{

  if(length(fastq_list)<3){
  filt= filterAndTrim(fwd= fastq_list$fastq,
                filt = filt_list,
                truncLen =cutting_param,
                multithread =cores,
                verbose=T,
                maxEE = maxEE,
                trimLeft = trimleft, ...)
  }
  if(length(fastq_list)>2){
   filt= filterAndTrim(fwd= fastq_list$fastq_fwd,
                  rev =fastq_list$fastq_rv,
                  filt = filt_list$fwd,
                  filt.rev = filt_list$rev,
                  truncLen = cutting_param,
                  multithread =cores,
                  verbose=T,
                  maxEE = maxEE,
                  trimLeft = trimleft, ...)
  }
  return(filt)
}
