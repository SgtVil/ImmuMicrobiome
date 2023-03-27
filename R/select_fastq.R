#' Select the fastq in the directory you want.
#'
#' @param names_list A list of names
#' @param fastq_list A list of fastq
#'
#' @return A list.
#' @export
#'
#' @examples No example
select_fastq = function(names_list, fastq_list){
  fastq_list$fastq_fwd= fastq_list$fastq_fwd[fastq_list$names%in% names_list]
  fastq_list$fastq_rv =  fastq_list$fastq_rv[fastq_list$names%in% names_list]
  fastq_list$names= fastq_list$names[fastq_list$names%in% names_list]
  return(fastq_list)
}
