#' Merge reads and remove chimera
#'
#' @param dd Object returned by \link{dereplicate} or \link{dada}.
#' @param filt_fastq_list The list of filtered fastq
#' @param cores Number of cores. Default=1
#' @param ... Parameter passed to \link{removeBimeraDenovo}
#'
#' @return
#' Return a count matrix.
#' @export
#'
#' @examples
#' No example.
merge_and_remove_chim= function(dd, filt_fastq_list, cores=1, ...){
  if(!is.list(filt_fastq_list)){
     pseudo_mergers <- dd %>%
     makeSequenceTable(.)%>%
     removeBimeraDenovo(., multithread=cores, method = "pooled", verbose = T)
 }
  if(is.list(filt_fastq_list)){
    message("merging an removing chimera from pair end")
    ddFs= dd$ddFs
    ddRs= dd$ddRs
    fwd= filt_fastq_list$fwd
    rv = filt_fastq_list$rev
    mergers <- mergePairs( dadaF = ddFs, derepF =fwd , dadaR =  ddRs,derepR = rv
                          , verbose = T, ...)

    pseudo_mergers <- mergers %>%
      makeSequenceTable(.)%>%
      removeBimeraDenovo(., multithread=cores, method = "pooled", verbose = T)
  }
  return(pseudo_mergers)
}
