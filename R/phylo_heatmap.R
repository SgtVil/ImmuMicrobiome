#'A function to plot heatmap nicely, using a phylogenetic distance.
#'@param physeq (Required). A phyloseq object.
#'@param top The top number of ASV you want to represent. Default is  \code{30}. It's usually the best number to plot.
#'@param distance Phylogenetic dinstance to use as dendrogramm. See \code{\link{beta_diversity}}.
#'@param taxa_rank The taxa rank at which you want to agregate abundances. Default is \code{"Genus"}.
#'@param split (Careful !). Numeric number by which you want the dendrogramm to be split.
#'Work only with ComplexHeatmap version 2.2.0.
#'@param factor_to_plot (Optional) Factor for the title.
#'@param color_vector A character string for color plottinh. Default is \code{c("white","#88CCAA","#771122")}.
#'@param taxa_names_par Graphical parameters. Default is \code{gpar(fontsize= 15, fontface= "bold.italic")}
#'@param label_names_par Graphical paramaters. Default is \code{gpar(fontsize=15, fontface="bold")}
#'@import ComplexHeatmap

#'@examples
#'library(phyloseq)
#'data(enterotype)
#'phylo_heatmap(physeq=enterotype, top=30, labels=NULL, distance="bray",
#'taxa_rank="Family")
#'@export
phylo_heatmap= function(physeq= physeq, labels= NULL, top= 30, factor_to_plot=NULL, distance="wunifrac", taxa_rank= "Genus",
                       split= NULL, color_vector= c("white","#88CCAA","#771122"),
                       taxa_names_par= gpar(fontsize= 15, fontface= "bold.italic"),
                       label_names_par= gpar(fontsize=15, fontface="bold")){

  if(is.null(access(physeq, "tax_table")) |is.null(access(physeq, "otu_table"))){
    stop("Need a phyloseq object")
  }
  if(is.null(access(physeq, "phy_tree"))){
    distance= "bray"
    print("No phylogenetic tree in this phyloseq object, bray-curtis distance selected.")
  }
  if(taxa_rank=="ASV"){
    H<- physeq%>%
      prune_taxa(names(sort(taxa_sums(.),TRUE)[1:top]), .)
    dend= hclust(phyloseq::distance(physeq, method=distance), "ward.D2")
  }
  H<- tax_glom(physeq, taxrank = taxa_rank)%>%
    prune_taxa(names(sort(taxa_sums(.),TRUE)[1:top]), .)
  dend= hclust(phyloseq::distance(physeq, method=distance), "ward.D2")
  if(taxa_are_rows(physeq)==TRUE){
    x= as.matrix(H@otu_table@.Data)
    print(" not reversed")
  }
  else{
    x= as.matrix(t(H@otu_table@.Data))
    print("reversed")
  }

  if(is.null(factor_to_plot)== TRUE & !is.null(labels)){

    # print("no factor")
    name= HeatmapAnnotation(Name= anno_text(as.matrix(sample_data(H)[, labels]), gp = gpar(fontsize=9)))
    h= Heatmap(as.matrix(scale(x)),
               cluster_rows = F, cluster_columns = as.dendrogram(dend), bottom_annotation = name,
               row_labels = tax_table(H)[, taxa_rank],
               row_names_gp = gpar(fontsize=15, fontface="bold.italic"),show_column_names = F,
               col = color_vector, column_dend_height = unit(3, "cm"), name = "Abundance",
               column_split = split, column_gap = unit(5, "mm"), column_title = NULL)
    print(h)
  } else if(is.null(labels) & !is.null(factor_to_plot)){
    # print("no labels")
    group= HeatmapAnnotation(factor_to_plot= anno_simple(as.matrix(sample_data(H)[,factor_to_plot])), show_legend = T)
    h= Heatmap(as.matrix(scale(x)),
               cluster_rows = F, cluster_columns = as.dendrogram(dend), top_annotation =  group,
               row_labels = H@tax_table[,taxa_rank],
               row_names_gp = gpar(taxa_names_par),show_column_names = F,
               col = color_vector, column_dend_height = unit(3, "cm"), name = "Abundance",
               column_split = split, column_gap = unit(5, "mm"), column_title = NULL)
    print(h)

  } else if(is.null(labels) & is.null(factor_to_plot)){
    # print("nothing")
    h= Heatmap(as.matrix(scale(x)),
               cluster_rows = F, cluster_columns = as.dendrogram(dend),
               row_labels = tax_table(H)[,taxa_rank],
               row_names_gp = gpar(fontsize=15, fontface="bold.italic"),show_column_names = F,
               col = color_vector, column_dend_height = unit(3, "cm"), name = "Abundance",
               column_split = split, column_gap = unit(5, "mm"), column_title = NULL)
    print(h)
  } else  if(!is.null(labels) & !is.null(factor_to_plot)){
    # print("both to plot")
    # print(dim(x))
    group= HeatmapAnnotation(factor_to_plot= as.matrix(sample_data(H)[,factor_to_plot]), show_annotation_name = T, name= factor_to_plot)
    name= HeatmapAnnotation(Name= anno_text(as.matrix(sample_data(H)[, labels]), gp = gpar(fontsize=9)))
    h= Heatmap(as.matrix(scale(x)),
               cluster_rows = F, cluster_columns = as.dendrogram(dend),
               top_annotation =  group, bottom_annotation = name,
               row_labels = H@tax_table[,taxa_rank],
               row_names_gp = gpar(taxa_names_par),show_column_names = F,
               col = color_vector, column_dend_height = unit(3, "cm"), name = "Abundance",
               column_split = split, column_gap = unit(5, "mm"), column_title = NULL,
               heatmap_legend_param = list(
                 group= list(title= factor_to_plot))
               )

    draw(h)

  }
}
