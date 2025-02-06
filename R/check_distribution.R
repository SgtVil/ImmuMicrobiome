#' Title
#'
#' @param df Dataframe to be used
#' @param values Values (numeric) to be used
#' @param factor Factor use for plotting
#' @param normalisation Normalisation method. So far the function accepts "log", "clr", "scaling", "sum scaling", "sqrt", "pareto", "crt", "vast", "level", "median centering", "gm centering" and "all".
#'"all" will give a facetted plot with all method depicted.
#'
#' @returns Returns a plot with the given normalisation method
#' @export
#'
#' @examples TBD
check_distribution = function(df, values, factor, normalisation = "none"){
source("~/Documents/ImmuMicrobiome/R/pareto_scale.R")

  norm_met= function(values, normalisation){
    switch(normalisation,
                   "log" = log(values),
                   "clr" = compositions::clr(values),
                   "scaling"= scale(values),
                   "sum scaling" = sum_scaling(values),
                   "sqrt" = sqrt(values),
                   "pareto"= pareto_scale(values),
                   "crt" = cubic_root(values),
                   "vast" = vast_scaling(values),
                   # "autoscale" = auto_scale,
                   "level"= level_scaling(values),
                   # "powerTransfrom" = car::powerTransform(values),
                   "median centering" = median_centering(values),
                   "gm centering"= geometric_mean(values))
  }

if(normalisation !="none" & normalisation != "all")  {
  # print("one")

  p=  df %>%
    group_by(across(all_of(factor))) %>%
    mutate(scaled_val= norm_met(eval(sym(values)), normalisation))  %>%
    ggplot(aes(scaled_val))+
    geom_density(aes_string(color=factor))+
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank())

  return(p)
}
  if( normalisation =="all"){

    print("all")
    p_list = list()

    for(normalisation in c("log", "clr", "scaling", "sum scaling", "sqrt", "pareto", "crt",
               "vast",  "level",  "median centering", "gm centering")){
# "powerTransform",
# print(normalisation)
      p_list[[normalisation]]=  df %>%
        group_by(across(all_of(factor))) %>%
        mutate(scaled_val= norm_met(eval(sym(values)), normalisation))  %>%
        ggplot(aes(scaled_val))+
        geom_density(aes_string(color=factor))+
        ggtitle(capitalize(normalisation))+
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5),
              axis.title.x = element_blank())

      # return(p_list)
    }

    p = gridExtra::grid.arrange(grobs=p_list )

    return(p)
  }

  else {
  # print("none")
  p=  df %>%
    ggplot(aes_string(values))+
    geom_density(aes_string(color=factor))
}

return(p)
}
