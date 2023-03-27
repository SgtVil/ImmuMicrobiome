#' Preprocess relative abundance.
#'
#' Process zero values and produce ratios for [slide_z()]. This function requires a \code{zero_treatment} value,
#' make sure to use [neg_dispersion] to see which method is the most convenient and accurate for your data.
#'
#' @param df A dataframe given by [seq_table()].
#' @param positive_sorted_sample Positive sorted sample pattern.
#' @param negative_sorted_sample First negative sorted sample pattern.
#' @param second_negative_sample Second negative sorted sample pattern.
#' @param zero_treatment Zero treatment to apply, see [random_generation()], [no_zero()] or [minimum_count()].
#'
#' @return Return a dataframe for each sample containing :
#'  \item{log10_abundance}{The \eqn{log10(pos*neg1)} ratio}
#'  \item{log2_ratio}{The \eqn{log2(pos1/neg1)}}
#'  \item{log2_neg_ratio} {The \eqn{log2(neg1/neg2)}}
#'  \item{log10_neg_abundance}{The \eqn{log10(neg1*neg2)}}
#'  \item{taxa} {A collapsed taxonomy}
#'
#' @export
#'
#' @examples No example
log_ratio <- function(df, positive_sorted_sample, negative_sorted_sample, second_negative_sample,
                      zero_treatment){
  tryCatch({
  pos = df[, which(grepl(colnames(df), pattern = positive_sorted_sample))]
  neg1 = df[, which(grepl(colnames(df), pattern = negative_sorted_sample))]
  neg2 = df[, which(grepl(colnames(df), pattern = second_negative_sample))]
# pos= df[,"pos"]
# neg1= df[,"neg1"]
# neg2= df[,"neg2"]
  # if(!is.null(dim(pos | neg1 | neg2))) { }

  if(zero_treatment == "random generation"){
    vec = pos!=0 & neg1!=0 & neg2!=0
    s = sum(c(length(pos[pos==0]),length(neg1[neg1==0]), length(neg2[neg2==0])))
    random= sapply(df[,2:4],function(x)quantile(x,seq(0, 1, 1/20))) [2,] %>%
      mean()
    # %>%
    #   round
    min= min(c(pos[pos!=0], neg1[neg1!=0], neg2[neg2!=0]))
    # random= ifelse(random)
    pos[pos==0] <- runif(length(pos[pos==0]),
                         min= min,
                         max= min +random)
    # %>% round
    neg1[neg1==0] <- runif(length(neg1[neg1==0]),
                           min= min,
                           max= min+random)
    # %>% round
    neg2[neg2==0] <- runif(length(neg2[neg2==0]),
                           min= min,
                           max= min+random)
    # %>% round
    message(paste(s, "asv equal to 0 were replaced by a random generation with a max of ", random))
  }

  if(zero_treatment =="min count"){
    pos[pos==0] <- min(pos[pos!=0])
    neg1[neg1==0] <- min(neg1[neg1!=0])
    neg2[neg2==0] <- min(neg2[neg2!=0])
  }
  if(zero_treatment=="no zero"){
    vec = pos!=0 & neg1!=0 & neg2!=0
    pos <- pos[vec]
    neg1 <- neg1[vec]
    neg2 <- neg2[vec]
    df = df[vec,]
  }

  log10_abundance <- log10(pos*neg1)
  log2_ratio <- log2(pos/neg1)
  log2_neg_ratio <- log2(neg1/neg2)
  log10_neg_abundance <- log10(neg1*neg2)
  taxa= as.matrix(rownames(df)) %>% set_colnames("taxa")
  tmp = as.data.frame(cbind( pos, neg1, neg2, log10_abundance, log2_ratio, log10_neg_abundance, log2_neg_ratio))

  df %>% dplyr::select(taxonomy, sample_id, new)%>% bind_cols(tmp) %>% bind_cols(taxa)%>%
    arrange(desc(log10_abundance))
}, finally = function(e)warning(paste("There is duplicated samples or missing fraction in:", colnames(df))))
}
