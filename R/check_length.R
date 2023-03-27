#' Check the length of reads in fastq.
#'
#' @param fastq To do
#' @param names TBD
#' @param cores TBD
#0'
#' @return TBD
#' @export
#'
#' @examples No example
check_length = function(fastq, names, cores=1){
 f= parallel::mclapply(fastq, seqTools::fastqq, mc.cores = cores)
 s = lapply(f, maxSeqLen)
 # min = lapply(f, minSeqLen)
 # colnames(s) = fastq[[3]]
 s = do.call("rbind", s)
 s =data.frame(len=s, names= names)
 s
}

