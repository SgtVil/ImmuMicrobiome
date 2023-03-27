search_mismatch = function(fastq_list_fwd, fastq_list_rev){
  fF <- FastqStreamer(fastq_list_fwd, n = 1e+06)
  fR <- FastqStreamer(fastq_list_rev, n = 1e+06)
  
  fqF <- suppressMessages(yield(fF, qualityType = "Auto"))
  fqR <- suppressMessages(yield(fR, qualityType = "Auto"))
  
  if(length(fqF) != length(fqR)){
    message(paste(fastq_list_fwd[[1]], "doesn't match rev"))
    return(fastq_list_fwd[[1]])
  }
}

