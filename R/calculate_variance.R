#' Calculate the variance for each factors.
#'
#' @param meta Metadata containing the factors to be used.
#' @param feat Features matrix
#'
#' @return
#' @import reshape2
#' @import rstatix
#' @export
#'
#' @examples
calculate_variance = function(meta, feat, cores=1){
  rownames(meta)= paste("sample", 1:length(meta[,1]))
  rownames(feat)= paste("sample", 1:length(feat[,1]))
  df = feat %>%
    as.data.frame %>%
    rownames_to_column()%>%
    inner_join(meta %>%
                 rownames_to_column(), by='rowname')%>%
    pivot_longer(cols = 2:(dim(feat)[2]+1), names_to = "features", values_to = "value")

  # var.tot = apply(feat, 2, function(x)var(x)) %>% as.data.frame
  mean.feat = apply(feat, 2, mean) %>%
    as.data.frame %>%
    rownames_to_column()%>%
    set_colnames(c("features", "mean.feat"))
  # x = feat[, x]
  var.grp = lapply(colnames(meta), function(i){
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
  names(var.grp)= colnames(meta)

  p.value = parallel:: mclapply(colnames(meta), function(i){
  if(length(unique(meta[,i]))>2){
   tmp= df  %>%
      drop_na(i)%>%
      group_by(features) %>%
      select_at(vars(features, value, i))%>%
      rstatix::kruskal_test(formula= as.formula(paste("value~", i))) %>%
      rstatix::adjust_pvalue(method = "fdr")%>%
      select(features, p, p.adj)%>%
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
      select(features, p, p.adj)%>%
      mutate(factor= i)%>%
      labelled::remove_attributes( 'all')
      # attributes(tmp)= NULL
  }
return(tmp)
  },mc.cores = cores)

  names(p.value)= colnames(meta)
  var.exp= reshape2::melt(var.grp) %>%
    # as.data.frame()%>%
    set_colnames(c("features", "variable", "var.grp", "factor")) %>%
    pivot_wider(names_from = factor, values_from = var.grp)
  # colnames(var.exp)= colnames(var.grp)

  p.value= do.call("rbind", p.value)

  return(list(variance= var.exp,
              p.value= p.value,
              mean.feat = mean.feat,
              all= var.grp))
}


