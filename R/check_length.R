#' Check the length of reads in fastq.
#'
#' @param fastq List of fastq
#' @param names Names of the samples
#' @param cores Number of threads to use. Default=1
#0'
#' @return TBD
#' @export
#' @return A dataframe with the length of sequences for each fastq.
#' @examples No example
check_length = function(fastq, names, cores=1){
 f= parallel::mclapply(fastq, seqTools::fastqq, mc.cores = cores)
 s = lapply(f, maxSeqLen)
 # min = lapply(f, minSeqLen)
 # colnames(s) = fastq[[3]]
 s = do.call("rbind", s)
 s = data.frame(len=s, names= names)
 s
}

