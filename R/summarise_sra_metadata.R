summarise_sra_metadata = function(dir, regex_match, host="human|mouse"){
  f= list.files(dir, full.names = T)
  files = map(f, read.csv) # Read files
  names(files)= sapply(strsplit(basename(f),"_"), "[", 1)

  # get column names
  column_names= map(files, colnames)
  com_colnames = Reduce(intersect, column_names)
  not_com_colnames = Reduce(union, column_names)
  not_com_colnames = not_com_colnames[not_com_colnames %notin% com_colnames]

  # get the mean bytes
  size= map_df(files, function(x){round(mean(x[,"Bytes"]/10^6), digits = 1)})

  # number of samples
  sample_nb = map_df(files, function(x)dim(x)[1])

  # Average length
  length = map_df(files, function(x)round(median(x[,"AvgSpotLen"]), digits = 1))

  # Instrument used
  instrument = sapply(files, "[", "Instrument")
  instrument = sapply(instrument, unique)

    if(length(instrument)!= length(files)) stop("not single instrument, need to revise this function")

  # is it human samples, humand and other, or no informations
  tmp=list()
  for(i in 1:length(files)){
    suppressWarnings({
   logic= grep(pattern = host, files[[i]], ignore.case = T)
    tmp[[i]] = files[[i]][ ,logic] %>%
                 str_extract_all(pattern = regex("human|mouse"), simplify = F)
    tmp[[i]] = unlist(tmp[[i]]) %>%
                  unique

    })
  }

  human_or_other=  map(tmp, function(x) case_when(
    length(x)==1 ~ "human",
    length(x)>1 ~ "human and other",
    length(x)==0 ~ "no match"
                      )) %>%
                          unlist
  # the name of the age variable
  tmp=list()
  for(i in 1:length(files)){
    suppressWarnings({
      logic= grep(pattern = "age", colnames(files[[i]]),ignore.case = T )
      tmp[[i]] = colnames(files[[i]])[logic]
         if(length(tmp[[i]])>=1) tmp[[i]]= paste(tmp[[i]], collapse = "|-|")
             if(length(tmp[[i]])==0) tmp[[i]]= "no data on age"
  })
  age= unlist(tmp)
  }
  # regex matching
  bool= map(files, grepl, pattern= regex_match)
      if(any(map(bool, sum)!=0)){
        message("The SRA metadata contains one of the regex_match, you should run clean_sra_data")
        for(i in 1:length(bool)){
          if(sum(bool[[i]])!=0) bool[[i]] = "regex found in this metadata"
          else bool[[i]] = "no regex found in this metadata"
        }
      }

  # export
  suppressWarnings(
    {df= as(rbind(size, sample_nb, length, human_or_other, instrument, age, regex_match= unlist(bool) ), "data.frame")
     rownames(df)= c("Mean size (Mo)", "Number of sample", "Sequencing length", "Type of samples", "Instrument used", "Name(s) of age variable", "Regex found or not")
      }
  )

  return(list(sum_metadata_sra= df, com_colnames= com_colnames, not_com_colnames= not_com_colnames))
}
