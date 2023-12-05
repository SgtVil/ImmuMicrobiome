#' Calculate fold change and p.value for your data.
#'
#' Accepts a metadata and a feature/microbiome matrix.
#' @param mat A features/ASV table with the sample in rows.
#' @param clinical_data A metadata containing the factors you want to test.
#' @param cores Number of cores you want to use. Default = 1.
#' @param ... Argument passed to \link{kruskal_test} or  \link{wilcoxon_test} depending on the number of groups to compare.
#'
#' @return
#' This function will return a list of data.frames. Each factor present in your `clinical_data` will be processed
#' and give an independent data.frame. Each dataframe contains:
#' - p = p value calculated by either Kruskal test or Wilcoxon test
#' - p.adj = p value corrected with FDR
#' - fold = fold change calculated as : log2(group1/group2)
#' @export
#'
#' @examples
#' fold= fold= fold_change(metabolomic[,1:4], features = metabolomic[,5:217], cores = 5)
#' head(fold$birth_type, 5)
fold_change = function(mat, clinical_data, cores = 2, ...){
  features = mat[,-clinical_data]
  clinical_data = mat[,clinical_data]

  rownames(clinical_data)= paste("sample", 1:length(clinical_data[,1]))
  rownames(features)= paste("sample", 1:length(features[,1]))

  df = features %>%
    as.data.frame %>%
    rownames_to_column()%>%
    inner_join(clinical_data %>%
                 rownames_to_column(), by='rowname')%>%
    pivot_longer(cols = 2:(dim(features)[2]+1), names_to = "features", values_to = "value")%>%
    as.data.frame



  p.value = parallel::mclapply(colnames(clinical_data), function(i){

      tmp= df  %>%
        drop_na(i)%>%
        group_by(features) %>%
        select_at(vars(features, value, i))%>%
        rstatix::wilcox_test(formula= as.formula(paste("value~", i)), p.adjust.method = "none", ...) %>%
        rstatix::adjust_pvalue(method = "fdr")%>%
        dplyr::select(p, p.adj)%>%
        mutate(factor= i)%>%
        labelled::remove_attributes( 'all')
      # attributes(tmp)= NULL

    return(tmp)
  }, mc.cores= cores)
  names(p.value)= colnames(clinical_data)

  fold.change= parallel::mclapply(colnames(clinical_data), function(i){
    if(length(unique(df[,i]))==2){
      comp=  df %>%
        drop_na(i)%>%
        # pull(!!sym(i))%>%
        distinct(!!sym(i)) %>%
        pull %>% sort

      tmp= df %>%
        drop_na(i)%>%
        group_by(features, !!sym(i)) %>%
        dplyr::select(features, !!sym(i), value)%>%
        summarise(mean_val=mean(value))%>%
        summarise(fold= log2(mean_val[2]/mean_val[1]),
                  mean_val= mean(mean_val))%>%
        ungroup%>%
        dplyr::select(features, fold) %>%
        mutate(comp = paste0(i, " : ","log2(", comp[2], "/", comp[1], ')'))

    } else if(length(unique(df[,i]))>2){
      comp=  df %>%
        drop_na(i)%>%
        # pull(!!sym(i))%>%
        distinct(!!sym(i)) %>%
        pull %>% sort

      tmp= df %>%
        drop_na(i)%>%
        group_by(features, !!sym(i)) %>%
        dplyr::select(features, !!sym(i), value)%>%
        summarise(mean_val=mean(value))%>%
        summarise(fold.1= log2(mean_val[2]/mean_val[1]),
                  fold.2= log2(mean_val[3]/mean_val[1]),
                  fold.3= log2(mean_val[3]/mean_val[2]))%>%
        ungroup%>%
        magrittr::set_colnames(value = c("features",  paste0(i, " : ","log2(", comp[2], "/", comp[1], ')'),
                               paste0(i, " : ","log2(", comp[3], "/", comp[1], ')'),
                               paste0(i, " : ","log2(", comp[3], "/", comp[2], ')')))%>%
        pivot_longer(!features, names_to = "comp", values_to = "fold" )
    }})

    mean.feat =  lapply(colnames(clinical_data), function(i){
      mean.feat = df %>%
      drop_na(i)%>%
      group_by(features) %>%
      select_at(vars(features, value, i))%>%
      summarise(mean_val = mean(value),
                median_val= median(value))%>%
        dplyr::select(median_val, mean_val)
      })



    fold.change = Map(cbind, mean.feat, p.value, fold.change)
    names(fold.change)= colnames(clinical_data)

  return(fold.change)
}
