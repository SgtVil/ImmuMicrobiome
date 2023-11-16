#' Plot the variance for two factors.
#'
#' @param variance An object returned by \link{calculate_variance}.
#' @param x Factor to be plotted as x.
#' @param y Factor to be plotted as y.
#' @param corrected Use p value adjusted by FDR or not. Default=F.
#' @param col Color vector.
#'
#' @return
#' A ggplot.
#'
#' @export
#'
#' @examples
#'
#' var = metabolomic %>%
#' dplyr::select(!child_id) %>%
#' calculate_variance(clinical_data = 1:3, cores = 1)
#' plot_xy_variance(var, x = "birth_type", y = "breastfeeding", corrected = F)
#'
plot_xy_variance = function(variance, x, y, corrected=F, col_vector= c("purple", "orange", "brown","grey"), max.overlap=20, size=3){

  # OLD WAY
  # ifelse( variance$p.value %>%
  #           filter(factor==fac)%>%
  #           pull(p.adj)<0.05, col[1], ifelse( variance$p.value %>%
  #                                               filter(factor==fac)%>%
  #                                               pull(p)<0.05, col[2], col[3]))
  if(corrected==F){
    col.val = variance$p.value %>%
      filter(factor==x | factor==y) %>%
      select(!p.adj) %>%
      pivot_wider(names_from = factor, values_from = p)

    tmp= ifelse(col.val[,x]<0.05 & col.val[,y]<0.05, "Both",
                no = ifelse(col.val[,x]<0.05 & col.val[,y]>0.05, colnames(col.val[,x]),
                            no = ifelse(col.val[,x]>0.05 & col.val[,y]<0.05, colnames(col.val[,y]),
                                        no="None")))
    labels = variance$p.value %>%
      # filter(factor==fac & p.adj<0.05)
      filter(p<0.05)
  } else {
    col.val = variance$p.value %>%
      filter(factor==x | factor==y) %>%
      select(!p) %>%
      pivot_wider(names_from = factor, values_from = p.adj)
    tmp= ifelse(col.val[,x]<0.05 & col.val[,y]<0.05, "Both",
                no = ifelse(col.val[,x]<0.05 & col.val[,y]>0.05, colnames(col.val[,x]),
                            no = ifelse(col.val[,x]>0.05 & col.val[,y]<0.05, colnames(col.val[,y]),
                                        no="None")))
    labels = variance$p.value %>%
      # filter(factor==fac & p.adj<0.05)
      filter(p.adj<0.05)
  }

col.val$col= tmp



  max.var= variance$variance %>%
    filter(variable=="var.exp")%>%
    summarise(across(c(!!sym(x), !!sym(y)), .fns=max))%>%
    max
  min.var= variance$variance %>%
    filter(variable=="var.exp")%>%
    summarise(across(c(!!sym(x), !!sym(y)), .fns=min))%>%
    min

  p=variance$variance %>%
    filter(variable=="var.exp")%>%
    left_join(col.val %>%
                select(features, col), by="features")%>%
    ggplot()+
    geom_point(aes(fill=col, !!sym(x), !!sym(y), size= variance$mean.feat$mean.feat), shape=21, alpha=0.7)+
    scale_fill_manual(values= col_vector,
                      name= "P values",
                      labels = c("Both factors", colnames(col.val[,x]), colnames(col.val[,y]), "None"),
                      breaks = c("Both", colnames(col.val[,x]), colnames(col.val[,y]), "None"))+
    # scale_fill_identity(name = "P values",
    #                     breaks = c("Both factors" = col_vector[1], x=col_vector[2], y= col_vector[2], "None"= col_vector[4]),
    #                     labels = c("Both factors", x, y, "None"),
    #                     guide = "legend")+
    geom_abline(slope = 1, intercept = 0, linetype="dashed")+
    ggrepel::geom_text_repel(data= variance$variance %>%
                       filter(features %in% labels$features & variable=="var.exp"),
                     aes(label= features, !!sym(x), !!sym(y)), max.overlaps = max.overlap, size=size, fontface="bold")+
    guides(size= "none", fill= guide_legend(override.aes = list(size=3)))+
    scale_y_continuous(limits = c(min.var, max.var), expand=c(0,0.005), labels=scales::percent, trans = "sqrt")+
    scale_x_continuous(limits = c(min.var, max.var), expand=c(0,0.005), labels=scales::percent, trans = "sqrt")
p
  return(p)

}

