#' Randomly select fastq files
#'
#' @param fastq_list List of fastq returned by \code{\link[vegan]{list_fastq}}
#' @param n Number of fastq to select
#'
#' @return
#' A list of fastq randomly selected.
#' @export
#'
#' @examples
#' TBD
rng_fastq= function(fastq_list, n=10){
  if(length(fastq_list$names)<n){
    ret = fastq_list
  } else {
    if(length(fastq_list)>2){
      rng= sample(1:length(fastq_list$names), n)
      ret = list(fastq_fwd= fastq_list$fastq_fwd[rng], fastq_rv= fastq_list$fastq_rv[rng], names= fastq_list$names[rng])
    } else {
      rng= sample(1:length(fastq_list$names), n)
      ret = list(fastq_fwd= fastq_list$fastq[rng], names= fastq_list$names[rng])
    }
  }
  return(ret)
}

