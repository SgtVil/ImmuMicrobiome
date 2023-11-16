#' Plot fold change or abundance as a ternary plot
#' @description
#' A wrapper function using \link{ggtern} function to produce a nice ternary plot.
#'
#'
#' @param fold A dataframe form the object returned by \link{fold_change} function.
#' Be aware that you need to specify which dataframe you want or use a `lapply` like approach. See examples for more.
#' @param adjusted.p.val Should the function filter features base on adjusted p values or not. Default =T
#' @param measure Should the function use fold change measure or the mean measures. Default = "fold.change".
#' @param stat Not implemented yet.
#'
#' @return
#'
#'
#' @examples
ternary_plot = function(fold, adjusted.p.val = T, measure="fold.change", stat="median"){
  if(length(unique(fold$comp))!=3) stop("ternary_plot is meant for comparisons between three factors")
require(ggtern)
#create polygons
polygons = data.frame(a= c(0.5, 0, 0, 0.5, 0, 0.5, 0.5, 0, 1, 1, 0.5, 1),
                         b= c(0, 0, 1, 0.5, 1, 0.5, .5, 1, 0, 0, 0.5, 1),
                         c= c(0.5, 1, 1, 0.5, 0, 0, .5, 1, 0, 1, 0.5, 0),
                         label= c(rep("MeCT 3", 4), rep("MeCT 2", 4), rep("MeCT 1", 4)))
  if(adjusted.p.val== T){
    df=  fold %>%
      filter(p.adj < 0.05) %>%
      pivot_wider(names_from = comp, values_from = fold.change)

  } else {
      filter(p < 0.05) %>%
      pivot_wider(names_from = comp, values_from = fold.change)

  }

 p= df %>%
   ggtern(aes(`MeCT 1`, `MeCT 2`, `MeCT 3`))+
   geom_polygon(data= polygons, aes(a, b, c, fill= label), alpha=0.3, show.legend = F)+
   scale_fill_manual(values=mect_col)+
   geom_crosshair_tern(data= cross,
                       aes(x = `MeCT 1`, y = `MeCT 2`, z = `MeCT 3`), linetype=2, )+
   geom_point(aes(color=ifelse(p.adj<0.05, "FDR < 0.05",
                               ifelse(p<0.01,  "nFDR < 0.01", "nFDR < 0.05")), size=-log10(p)))+
   geom_text(data = df,  aes(label= features), size=3,  check_overlap = T, col="black",
             position = position_nudge_tern(0.01, .01,.01))+
   scale_color_manual(values=c("lightsalmon","#205062", "#42858C"),name="p values")+
   scale_R_continuous(labels = waiver())+
   Rlab(label = "MeCT 3")+
   Tlab("MeCT 2")+
   Llab("MeCT 1")+
   theme_bw()+
   theme_nomask()+
   theme(axis.title = element_text(size=15))
 return(p)
}
