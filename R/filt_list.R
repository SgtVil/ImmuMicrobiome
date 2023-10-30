#' Create a list of filtered fastq names.
#'
#' @param fastq_list A list of fastq returned by \link{list_fastq}
#'
#' @return
#' A list of filtered fastq.
#' @export
#'
#' @examples
#' No example.
filt_list <- function(fastq_list){
 path= unique(dirname(fastq_list[[1]]))
  if(length(fastq_list)<3){
    filt=  file.path(path, "filtered",
                paste0(fastq_list[["names"]] ,"_filt.fastq.gz"))
    return(filt)
    }
  if(length(fastq_list)>2){
  fwd= file.path(path, "filtered",
             paste0(fastq_list[["names"]] ,"_fwd_filt.fastq.gz"))
  rev= file.path(path, "filtered",
             paste0(fastq_list[["names"]] ,"_rv_filt.fastq.gz"))
  return(list(fwd=fwd, rev=rev))
 }

}
