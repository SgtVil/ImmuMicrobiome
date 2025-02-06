#' Title
#'
#' @param data A dataframe
#' @param formula The right side of the formula you want to apply to permanova
#' @param clinical_data Column names or number that are metadata
#' @param distance Distance to use, using vegan distances
#' @param by Either margin or terms
#'
#' @returns Returns a plot and the p-values of the permanova.
#' @export
#'
#' @examples
permutation_analysis = function(data, formula, clinical_data, distance="euclidean",  by="margin"){
`%notin%` = negate(`%in%`)
  # "as(otu_table(reverseASV(tmp)), 'matrix') ~ Constipation+Sex+Age_sampling"
  res = list()
  R2 = list()
  p.adonis = list()

  # formula= terms(from) %>% delete.response()
  form = as.formula(paste("data[,-clinical_data]", formula) )

    res = adonis2(form,
                       data = as(data, "data.frame"),
                       permutations = 999, na.action = na.exclude,
                       method = distance, by="margin")
    # res = adonis2(formula = flow_cyto[,3:116] ~ dm_intervention + dm_alter +recruitment_site, data = flow_cyto[, -2], by="margin")

    R2= data.frame(factor= rownames(res), R2= res$R2, p.value =  res$`Pr(>F)`)

   p= R2 %>%
      filter(factor %notin% c("Total", "Residual") ) %>%
      # drop_na() %>%
      # mutate(omic= factor(omic, levels=c("taxMG", "funMG", "taxMT", "funMT", "taxMP", "funMP","MM")),
      #        comp= gsub("Age_sampling","Age", comp),
      #        comp= factor(comp, levels=c("Sex", "Constipation", "Age")))%>%
      add_significance(p.col = "p.value")%>%
      ggplot(aes(factor, R2, fill=R2,  height=sqrt(1-p.value),
                 width=sqrt(1-p.value)))+
      geom_tile()+
      geom_text(aes(label= p.value.signif))+
      scale_fill_gradient(low="grey", high="darkslategray4")+
      # coord_flip()+
      theme(axis.title = element_blank(),
            text=element_text(size=15))

   return(list(p, R2))
}

#   R2.tax.cof=  R2
#   rownames(R2.tax.cof) = rownames(res[[1]])
#   colnames(R2.tax.cof)= c("taxMG", "taxMT")
#   res.adonis.tax.cof =p.adonis
#
#   colnames(res.adonis.tax.cof)= c("taxMG", "taxMT")
#
#   rownames(res.adonis.tax.cof) = rownames(res[[1]])
# }
