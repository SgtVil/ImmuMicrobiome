#' Find primer sequences in the fastq.
#'
#' @description
#' A short function to verify if the primers are still in the fastq. Usually primers are kept in the reads after the removal of adapters.
#'
#' @param fastq_list List of fastq
#' @param fwd_primer String of the forward primer
#' @param rv_primer String of the reverse primer
#' @param max.mismatch Max number of mismatches. This needs to be adapted depending on the number of wobble nucleotides from you primers
#' @param allow.indels Indels stands for insertion-deletion. Refers to \link{matchPattern} for more information on this parameter.
#' @param pair_end Logical.
#' @param fixed See \link{matchPattern} for full informations.
#'
#' @return
#' An object \link{vcountPattern}.
#' @export
#'
#' @examples
#' No example.
find_primer <- function(fastq_list, fwd_primer, rv_primer, max.mismatch, allow.indels=FALSE, pair_end= T, fixed=F){
  require(Biostrings)
  require(ShortRead)
  if(!is.character(fwd_primer)){
    stop(message("need a character vector surch as : ACGTYNNT"))
  }
  if(pair_end==TRUE){
    if (is.null(rv_primer)){
      rd= ShortRead::readFastq(fastq_list[[1]][1]) %>% sread
      f= vcountPattern(pattern = fwd_primer, subject = , max.mismatch = max.mismatch,
                      with.indels = allow.indels)
      return(f)
    } else{
      rd_fwd= ShortRead::readFastq(fastq_list[[1]][1])%>% sread
      rd_rv= ShortRead::readFastq(fastq_list[[2]][1]) %>% sread
      fwd= vcountPattern(pattern = fwd_primer, subject = rd_fwd, max.mismatch = max.mismatch,
                        with.indels = allow.indels, fixed = fixed)
      rv= vcountPattern(pattern = rv_primer, subject = rd_rv, max.mismatch = max.mismatch,
                       with.indels = allow.indels, fixed = fixed)
      return(list(fwd, rv))
    }
  }
  if(pair_end== FALSE){
    if (is.null(rv_primer)){
      rd= ShortRead::readFastq(fastq_list[[1]][1]) %>% sread
      f= vcountPattern(pattern = fwd_primer, subject = , max.mismatch = max.mismatch,
                      with.indels = allow.indels, fixed = fixed)
      return(f)
    } else{
      rd_fwd= ShortRead::readFastq(fastq_list[[1]][1])%>% sread
      rd_rv= ShortRead::readFastq(fastq_list[[1]][1]) %>% sread
      fwd= vcountPattern(pattern = fwd_primer, subject = rd_fwd, max.mismatch = max.mismatch,
                        with.indels = allow.indels, fixed = fixed)
      rv= vcountPattern(pattern = rv_primer, subject = rd_rv, max.mismatch = max.mismatch,
                       with.indels = allow.indels, fixed = fixed)
      return(list(fwd, rv))
    }
  }
}


