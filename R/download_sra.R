download_sra = function(accession_dir, output_dir, string="_SRR_Acc_List.txt", path_to_sra) {

  if(!dir.exists(output_dir)) dir.create(output_dir)
  f = list.files(accession_dir, pattern = string )
  if(!dir.exists("bash_files")) dir.create(paste0(getwd(), "/bash_files"))
  func_fetch<- paste0(path_to_sra, "/prefetch --option-file")
  for(i in f){
    accession<- paste0( accession_dir, "/", i )
    outdir= paste0(output_dir, "/", strsplit(i, split = "_")[[1]][1])
    cmd1<- paste(func_fetch, accession , "--output-directory",  outdir, "--progress") # here goes the function
     # cmd1 = paste0("cat ", accession_dir, "/", i)
    bash_out = dirname("/home/remy/Documents/Meta_analyse_ECOBIOTIC/Accession_list ")
    write(cmd1, file= paste0(bash_out,"/bash_files/", i, ".sh")) # write bash files
  }
  bash_cmd= paste0("bash /home/remy/test_bash.sh /home/remy/Documents/Meta_analyse_ECOBIOTIC/bash_files/*")
    system(bash_cmd)
}
