#' Title
#' Find differential abundance within a dataset.
#'
#' This is a wrapper function that use [ALDEx2] package to find the taxa that are differentially abundant
#' between two groups. You have to provide a phyloseq object with only two groups in order to use this function.
#'
#' @param physeq Phyloseq object to compute analysis.
#' @param group The group to perform analysis. Make sure the vector is of length 2. If not please subset your data.
#' @param paired Are the data paired.
#' @param col1 First column used to name the taxa in the plot. Default is Genus level.
#' @param col2 Second column used to name the taxa in the plot. Default is the last column of your
#' tax_table. It is recommended to use an ASV or OTU number for each taxa.
#' @return A list of three elements.
#' \item{all_features}{A dataframe containing all taxa tested with their respective values.}
#' \item{signif_features}{A dataframe containing only taxa significant for eBH test.}
#' \item{barplot}{Differential abundance plotted witha barplot}
#' \item{volcano}{Differential abundance plotted witha volcano plot}
#' @export
#' @import ggrepel
#' @importFrom ALDEx2 aldex.clr
#' @importFrom ALDEx2 aldex.effect
#' @importFrom ALDEx2 aldex.ttest
#'
differential_abundance <- function(physeq, group, paired= F, col1="Genus", col2=NULL, cores=F, plot=F){
  theme_set(theme_minimal()+
              theme(panel.grid = element_blank(), legend.position = "bottom",
              axis.title.y = element_blank(), axis.text.y = element_text(size = 20),
              legend.text = element_text(size = 20), axis.title.x = element_blank(),
              axis.text.x = element_text(size= 20),
              title = element_text(size = 20, face="bold")))
  # taxa_names(physeq)= paste(tax_table(physeq)[,"Genus"], 1:dim(tax_table(physeq))[1])

  if(taxa_are_rows(physeq)){
    otu= as.matrix(otu_table(physeq))
  } else {
    otu= as.matrix(t(otu_table(physeq)))
  }

  conds= as_vector(sample_data(physeq)[,group])
  x_clr= aldex.clr(otu, conds = conds)
  x= aldex.ttest(x_clr, paired.test = paired)
  x1 = aldex.effect(x_clr, useMC = cores )
  x_fin= data.frame(x, x1)

  x= x_fin[x_fin$we.eBH<0.05 | x_fin$wi.eBH<0.05,]

 name= unique(conds)

 if(plot==T){
   p1=  ggplot(x, aes(y= diff.btw,
                      x= reorder(rownames(x), X = diff.btw),
                      fill=diff.btw<0))+
     geom_bar(stat = "identity", )+
     coord_flip()+
     scale_fill_manual(labels= c(name[1], name[2]), name=NULL,
                       values = c(col1,col2))


   p2= x_fin %>%
     rownames_to_column("taxa")%>%
     # separate(taxa, into = c("Family", "Genus"), sep = "_", remove = F)%>%
     ggplot(aes(diff.btw, -log10(we.eBH)))+
     geom_point(aes(size=rab.all), fill=
                  ifelse(x_fin$diff.btw<0 & x_fin$we.eBH<0.05, col1,
                         ifelse(x_fin$diff.btw>0 & x_fin$we.eBH<0.05, col2,"black")), shape=21)+
     geom_hline(yintercept = -log10(0.05))+
     geom_label_repel(aes(label=ifelse(we.eBH<0.05, taxa, "")), max.overlaps = 40)+
     scale_fill_manual(labels= c(name[1], name[2]), name=NULL,
                       values = c(col1,col2))+
     labs(y= "-log10(p.value) BH corrected", size="Mean relative abundance")

   return(list(all_features= x_fin, signif_features= x, barplot=p1, volcano= p2))
 } else {
   return(list(all_features= x_fin, signif_features= x))
 }
 }



