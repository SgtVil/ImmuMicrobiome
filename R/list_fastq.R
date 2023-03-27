#' Make the list of fastq files for the rest of the proceedings.
#'
#' This will enable other wrapper functions for dada2 pipeline.
#'
#' @param path Path of selected samples
#' @param pattern Pattern for forward and reverse fastq if paired sample analysis, single pattern for single end.
#' @param separator Separator by which you want to separate the filenames to obtain sample names. Default = "_".
#' @param level Which level of the separated filenames you want to keep as sample names.
#'
#' @return A list of fwd, rv and names.
#' @export
#'
#' @examples No example
list_fastq= function(path, pattern, separator="_", level){

f=list.files(path)
# is= sum(grepl(pattern = pattern, f))
  if(length(pattern)==1){
   f= list.files(path, full.names = T)
   f= f[grepl(pattern = "fastq|fastq.gz", x = basename(f))]
   n= file.exists(f) %>% sum()
   message(paste(n, "fastq are single end or merged"))

   tmp <- str_split(basename(f), pattern = separator, simplify = T)[,level]

   return(list(fastq=f, names= tmp))
 } else {
   fwd=  list.files(path=path, pattern =  pattern[1], full.names = T)
   rv= list.files(path = path, pattern =  pattern[2],full.names = T)
   fwd= fwd[grepl(pattern = "fastq|fastq.gz", x = basename(fwd))]
   rv= rv[grepl(pattern = "fastq|fastq.gz", x = basename(rv))]

   n1= file.exists(fwd) %>% sum()
   n2= file.exists(rv) %>% sum()
   message(paste(n1, n2, "fastq are double end"))

   tmp1 <- str_split(basename(fwd), pattern=separator, simplify = T)[,level]
   tmp2 <- str_split(basename(rv), pattern=separator, simplify = T)[,level]

   return(list(fastq_fwd= fwd, fastq_rv= rv, names= tmp1))

 }

}
