#' Title Wrapper function to analyse alpha diversity.
#' @description
#' This function will rapidly plot alpha diversity in ggplot2.
#'
#' @param physeq A phyloseq object
#' @param measure The type of alpha diversity measure you wish to use. These are directly taken from \link{estimate_richness}.
#' @param x String character to display on X axis
#' @param group String character to color/fill by.
#' @param plot_type For the moment the function accepts "line" or "boxplot". Default = "boxplot"
#' @param stat Logical. Show stats or not. Default= FALSE.
#' @param check_depth Logical. Only for integers count, doesn't work with relative abundance. Will create two plots :
#' - a first one being a boxplot based ggplot with depth as size of the jitters
#' - a second one being a dotplot with alpha diversity measure on Y and the depth on X.
#' @param size Numeric value for the dot size and line width.
#'
#' @return Return a ggplot graphic.
#' @export
#'
#' @examples
#' data(enterotype)
#'
#' alpha_diversity(enterotype, measure="Shannon", x="Enterotype", group="SeqTech", plot_type="boxplot")
#' # same analysis but with check_depth
#' alpha_diversity(GlobalPatterns, measure="Shannon", x="SampleType", plot_type="boxplot", check_depth=T)


alpha_diversity= function(physeq, measure= "Shannon", x, group=NULL, plot_type="boxplot", color_vector=  c("brown","darkgreen","orange","violet"), stat=FALSE, check_depth=F, size=1.5 ){
  ad = estimate_richness(physeq, measures = measure)
  ad = cbind(as(sample_data(physeq), "data.frame"), ad)
  ad[, x] = factor(ad[,x])

  if(plot_type=="boxplot"){
    p = ggplot(ad,aes_string(y=measure, x))+
      geom_boxplot(alpha=0, size=1.5)+
      geom_jitter(aes_string(fill=group), position = position_jitterdodge(jitter.width = 0.25), shape=21, size=size)+
      scale_fill_manual(values = color_vector)
    if(stat==TRUE){
      p= p+
        stat_compare_means()
    }
  }

  if(plot_type== "line"){
    p = ggplot(ad,aes_string(y=measure, x))+
      geom_jitter(aes_string(color=group))+
      stat_summary(geom = "line", aes_string(color=group, group=group), linewidth=size)+
      stat_summary(geom="pointrange", aes_string(color=group, group=group), linewidth=size, size=size)+
      scale_fill_manual(values = color_vector)
    if(stat==TRUE){
      p= p+
        stat_compare_means(aes_string(group=group))
    }
  }
  if(check_depth==T){
    ad$depth = sample_sums(physeq)
   p1= ggplot(ad, aes_string(y=measure, x))+
      geom_boxplot(alpha=0, size=1.5)+
      geom_jitter(aes_string(fill=group, size= "depth"), position = position_jitterdodge(jitter.width = 0.25), shape=21)+
     scale_fill_manual(values = color_vector)
   p2= ggplot(ad, aes_string(y=measure, "depth", fill=group))+
     geom_point(shape=21)+
     scale_fill_manual(values = color_vector)
     # facet_wrap(facets = group)
   p= ggarrange(p1, p2, common.legend = T)
  }
  return(p)
}

