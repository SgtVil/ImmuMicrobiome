clean_sra_metadata= function(dir, regex_match, output_path="./curated_metadata"){
  f= list.files(dir, full.names = T)
  files = map(f, read.csv) # Read files
  names(files)= sapply(strsplit(basename(f),"_"), "[", 1)
  
  #mutate and filter any rows from any column that match our regex match, aka shotgun data
  bool= map(files, grepl, pattern= regex_match) %>%
    map(sum) %>% unlist
  tmp=files[[which(bool==1)]]

  new_data = tmp%>%
    rowwise %>% 
    mutate(regex_match = any(str_detect(c_across(is.character), regex(regex_match)), na.rm = FALSE)) %>% 
    filter(!regex_match)
  if(dim(new_data)[1]==0) message("No samples remaining after curation")
  dir.create(output_path)
  write.csv(new_data, file= paste0(output_path, "/", names(files)[which(bool==1)], "_SraRunTable.txt"))

 }
    
