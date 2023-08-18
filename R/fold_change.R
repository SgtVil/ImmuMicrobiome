#' Calculate fold change and p.value for your data.
#'
#' Accepts a metadata and a feature/microbiome matrix.
#' @param meta A metadata containing the factors you want to test.
#' @param feat A features/ASV table with the sample in rows.
#' @param cores Number of cores you want to use. Default = 1.
#'
#' @return
#' @export
#'
#' @examples
fold_change = function(meta, feat, cores = 2){

  rownames(meta)= paste("sample", 1:length(meta[,1]))
  rownames(feat)= paste("sample", 1:length(feat[,1]))
  df = feat %>%
    as.data.frame %>%
    rownames_to_column()%>%
    inner_join(meta %>%
                 rownames_to_column(), by='rowname')%>%
    pivot_longer(cols = 2:(dim(feat)[2]+1), names_to = "features", values_to = "value")

  p.value = parallel::mclapply(colnames(meta), function(i){
    if(length(unique(meta[,i]))>2){
      tmp= df  %>%
        drop_na(i)%>%
        group_by(features) %>%
        select_at(vars(features, value, i))%>%
        rstatix::kruskal_test(formula= as.formula(paste("value~", i))) %>%
        rstatix::adjust_pvalue(method = "fdr")%>%
        select(features,  p, p.adj)%>%
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
  }, mc.cores= cores)
  names(p.value)= colnames(meta)
  fold.change= parallel::mclapply(colnames(meta), function(i){
  comp=  df %>%
      drop_na(i)%>%
      # pull(!!sym(i))%>%
      distinct(!!sym(i)) %>%
      pull %>% sort
 df= df %>%
    drop_na(i)%>%
      group_by(features, !!sym(i)) %>%
      select(features, !!sym(i), value)%>%
      summarise(mean_val=mean(value))%>%
      summarise(fold= log2(mean_val[2]/mean_val[1]))%>%
      ungroup%>%
      select(fold) %>%
      mutate(comp = paste0(i, " : ","log2(", comp[2], "/", comp[1], ')'))

    return(df)


  })
  names(fold.change)= colnames(meta)

  tmp = Map(cbind, p.value, fold.change)

 return(tmp)
}
