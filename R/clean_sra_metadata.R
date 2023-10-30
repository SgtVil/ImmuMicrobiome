#' Function to remove sample from sra files.
#' @description
#' SRA repository can contain metagenomic files, metatranscriptomic or 16S files. Depending of which type of files you wan't to work
#' you can use a regex to remove unwanted samples.
#'
#'
#' @param dir Metadata directory
#' @param regex_match The regex needed to remove unwanted files.
#' @param output_path Directory for the new files. Default = "./curated_metadata".
#'
#' @return
#' This function will search for your regex and remove the rows matching your regex. It will create curated files in a new directory.
#' @export
#'
#' @examples
#' No example.
clean_sra_metadata= function(dir, regex_match, output_path="./curated_metadata"){
  f= list.files(dir, full.names = T)
  files = purrr::map(f, read.csv) # Read files
  names(files)= sapply(strsplit(basename(f),"_"), "[", 1)

  #mutate and filter any rows from any column that match our regex match, aka shotgun data
  bool=  purrr::map(files, grepl, pattern= regex_match) %>%
    purrr::map(sum) %>% unlist
  tmp=files[[which(bool==1)]]

  new_data = tmp%>%
    rowwise %>%
    mutate(regex_match = any(str_detect(c_across(is.character), regex(regex_match)), na.rm = FALSE)) %>%
    filter(!regex_match)
  if(dim(new_data)[1]==0) message("No samples remaining after curation")
  dir.create(output_path)
  write.csv(new_data, file= paste0(output_path, "/", names(files)[which(bool==1)], "_SraRunTable.txt"))

 }

