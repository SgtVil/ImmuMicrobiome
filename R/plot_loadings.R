#' Title
#'
#' @param ml_obj An object of class ml or a list of ml class
#'
#' @return A ggplot based graphic
#' @export
#'
#' @examples
plot_loadings= function(ml_obj, top=20){

  #prepare plotting function
  plot_fun1 = function(x)ggplot(x, aes(fct_reorder(rowname, cumulated_importance), x=loading))+
    geom_col()+
    facet_wrap(~group)+
    theme_classic()

  #prepare data for plotting
  if(is.list(ml_obj)){
    loads= lapply(ml_obj, slot, "bacterial_loading")
    loads= lapply(loads, function(x){
      x %>%
        mutate(cumulated_importance= apply(., 1, sum)) %>%
        arrange(desc(cumulated_importance)) %>%
        top_n(n=top)  %>%
        rownames_to_column()%>%
        pivot_longer(names_to = "group", values_to = "loading", cols = !rowname & !cumulated_importance)
    })
#plot for the list of ml_obj
    lapply(loads, function(x){
      if(anyNA(x) | dim(x)[1]==0){
        message("no loadings for this iteration")
      } else {
        plot_fun1(x)
      }
    })
  } else {
    loads = slot(ml_obj, "bacterial_loading")
    loads= loads %>%
      mutate(cumulated_importance= apply(., 1, sum)) %>%
      arrange(desc(cumulated_importance)) %>%
      top_n(top)  %>%
      rownames_to_column()%>%
      pivot_longer(names_to = "group", values_to = "loading", cols = !rowname & !cumulated_importance)

    plot_fun1(loads)
  }
}



