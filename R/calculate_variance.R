#' Calculate the variance for each factors.
#'
#' @description
#' Calculate the variance of a given clinical (discrete value) on a given feature. This function allows you to quantify the variance
#' explained by a given factor for each features with the formula: \deqn{ explained variance= 1 - group variance / total variance}
#' This function will take count of missing data. If you have any `NA` for a given factor, the total variance for each feature will be calcultaed without the `NAs`.
#' @param mat Features matrix
#' @param clinical_data Metadata containing the factors to be used.
#' @param cores Number of threads to use. Default = 1.
#'
#' @return
#' A list of data.frames :
#' - variance: data.frame containing the values of variance for each factors with the var.tot (total variance of a features),
#' the var.grp (variance explained by the groups in each factor) and var.exp (explained variance by the groups out of the total variance).
#' - p.value = p value calculated by either Kruskal test or Wilcoxon test and corrected with FDR
#' - mean.feat = Mean value for the features
#' - all = a data.frame with var.tot, mean.feat, var.grp and var.exp
#' @import reshape2
#' @importFrom rstatix wilcox_test
#' @importFrom rstatix t_test
#' @importFrom rstatix function
#' @export
#'
#' @examples
#' var = calculate_variance(metabolomic[,1:4], feature = metabolomic[,5:217])
#' # show variance data.frame
#'  head(var$variance, 5)
#'
#'  # show p.value data.frame
#'   head(var$variance, 5)
#'
#'   # show mean.feat data.frame
#'    head(var$variance, 5)
calculate_variance = function(mat, clinical_data, cores=1){
  features = mat[,-clinical_data]
  clinical_data = mat[,clinical_data]

  rownames(clinical_data)= paste("sample", 1:length(clinical_data[,1]))
  rownames(features)= paste("sample", 1:length(features[,1]))
  df = features %>%
    as.data.frame %>%
    rownames_to_column()%>%
    inner_join(clinical_data %>%
                 rownames_to_column(), by='rowname')%>%
    pivot_longer(cols = 2:(dim(features)[2]+1), names_to = "features", values_to = "value")

  # var.tot = apply(feat, 2, function(x)var(x)) %>% as.data.frame
  mean.feat = apply(features, 2, mean) %>%
    as.data.frame %>%
    rownames_to_column()%>%
    magrittr::set_colnames(c("features", "mean.feat"))
  # x = feat[, x]
  var.grp = lapply(colnames(clinical_data), function(i){
    df%>%
      drop_na(i)%>%
      group_by(features)%>%
      mutate(var.tot= var(value),
             mean.feat= mean(value))%>%
      group_by_at(vars(i, features, var.tot, mean.feat))%>%
      summarise(var.tmp= sum((value - mean(value))^2),
                n= n())%>%
      ungroup %>%
      group_by(features, var.tot, mean.feat) %>%
      summarise(var.grp = sum(var.tmp)/sum(n))%>%
      ungroup %>%
      mutate(var.exp= 1-var.grp/var.tot)
    # return(var.grp)
  })
  names(var.grp)= colnames(clinical_data)

  p.value = parallel:: mclapply(colnames(clinical_data), function(i){
  if(length(unique(clinical_data[,i]))>2){
    tmp= df  %>%
      drop_na(i)%>%
      group_by(features) %>%
      select_at(vars(features, value, i))%>%
      rstatix::kruskal_test(formula= as.formula(paste("value~", i))) %>%
      rstatix::adjust_pvalue(method = "fdr")%>%
      dplyr::select(features, p, p.adj)%>%
      mutate(factor= i)%>%
      labelled::remove_attributes('all')
   # attributes(tmp)= NULL
  } else {
    tmp= df  %>%
      drop_na(i)%>%
      group_by(features) %>%
      select_at(vars(features, value, i))%>%
      rstatix::wilcox_test(formula= as.formula(paste("value~", i))) %>%
      rstatix::adjust_pvalue(method = "fdr")%>%
      dplyr::select(features, p, p.adj)%>%
      mutate(factor= i)%>%
      labelled::remove_attributes( 'all')
      # attributes(tmp)= NULL
  }
return(tmp)
  }, mc.cores = cores)

  names(p.value)= colnames(clinical_data)
  var.exp= reshape2::melt(var.grp) %>%
    # as.data.frame()%>%
    magrittr::set_colnames(c("features", "variable", "var.grp", "factor")) %>%
    pivot_wider(names_from = factor, values_from = var.grp)


  p.value= do.call("rbind", p.value)
  # p =  p.value %>%
  #   dplyr::select(features, p, factor)%>%
  #   pivot_wider(names_from = factor, values_from = p)
  # p.adj =  p.value %>%
  #   dplyr::select(features, p, factor)%>%
  #   pivot_wider(names_from = factor, values_from = p)

  return(list(variance= var.exp,
              p.value= p.value,
              mean.feat = mean.feat,
              all= var.grp))
}


