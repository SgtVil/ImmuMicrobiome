#' IgASeq table
#'@description Prepare phyloseq object to perform the slide_z function.
#'
#'@param physeq A phyloseq object to analyse your data. (\code{\link{phyloseq-class}})
#'@param sample_name Character vector to identify the unique sample name from which you performed
#'the sortings
#'@param sorting_names Character vector containing the variable for the sorting fractions.
#'@param cols_to_keep The columns of \code{\link{phyloseq::sample_data-class}}. Default = "all", better to keep all of the columns.

#' @return A list equal to the length of unique samples that originated the sorted samples. The list contain the taxa collapsed to a single column (# as separator),
#' the sample_data() from the given phyloseq object and the pos, neg1 and neg2 relative abundance.
#'@export
seq_table <- function(physeq, sample_name, sorting_names, cols_to_keep="all") {
  sample_meta = sample_data(physeq) %>%
                 as.matrix() %>%
                  as.data.frame()
  sample_list = list()

  if(dim(otu_table(physeq))[1] < dim(otu_table(physeq))[2]){
    count_table = otu_table(physeq) %>%
      t %>%
      as.data.frame
    message("otu_table has been transposed")
  } else {
    count_table = otu_table(physeq) %>%
      as.data.frame()
  }

  if(is.null(access(physeq, "tax_table"))){
    stop("Need a phyloseq object")}
  if(sum(duplicated(sample_meta[, sorting_names]))>1) {
    warning("There is duplicated samples in your data, make sure to remove or change names of duplicates")}

  if (dim(sample_data(physeq))[1]<2){
    stop("no samples left after subseting to sorted samples")
  }
  if(cols_to_keep=="all"){
    new_cols= unite(as.data.frame(sample_meta), col = "all", sep="#")
  } else {
    new_cols= sample_meta %>%
      as.data.frame()%>%
      dplyr::select(all_of(cols_to_keep))%>%
      unite(col = "new",sep = "#")
  }

  for(i in unique(sample_meta[, sample_name])){

    iD <- rownames(sample_meta[sample_meta[, sample_name]==i,])


    tmp = count_table[,iD]

    if(is.null(dim(tmp))){
      cat("sample",i, "is alone", "\n")
      next
    }
    if(anyNA(tmp) | any(apply(tmp, 2, sum) == 0)){
      cat("One sample belonging to", i, "has no reads\n")
      next
    }

    if(dim(tmp)[2] > 2 & dim(tmp)[2]<5){
    tmp2 = data.frame(taxonomy = unite(col = taxonomy, as.data.frame(tax_table(physeq)), 1:dim(tax_table(physeq))[2], sep = "#"))
    tmp3 = cbind(tmp2, tmp)
    tmp3= tmp3[apply(tmp3[,-1],1, function(x)sum(x>0)>1)%>% which,]
    # %>%   set_colnames(c("taxonomy", "pos", "neg1", "neg2"))
    tmp3$sample_id= i
    tmp3$new= new_cols[rownames(new_cols) %in% iD[1],]
    sample_list[[i]]= tmp3
      }
    if(dim(tmp)[2] > 2 & dim(tmp)[2]>4){
      tmp2 = data.frame(taxonomy = unite(col = taxonomy, as.data.frame(tax_table(physeq)), 1:dim(tax_table(physeq))[2]))
      tmp3 = cbind(tmp2, tmp)
      tmp3= tmp3[apply(tmp3[,-1],1, function(x)sum(x>0)>1)%>% which,]
      # %>% set_colnames(c("taxonomy", "pos", "neg1", "neg2"))
      tmp3$sample_id= i
      tmp3$new= new_cols[rownames(new_cols) %in% iD[1],]
      sample_list[[i]]= tmp3
    }
    }
  return(sample_list)
}

