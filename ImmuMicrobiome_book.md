Introduction to the ImmuMicrobiome package
================
Rémy Villette
2023-10-26

# Introduction

<div style="text-align: justify">

This package is based on the code of Dr. Villette and Dr. Larsen. The
package is written and maintained by Dr. Villette. This package is meant
to facilitate microbiome exploration and ensuring nice plotting.  
This package covers :

- The dada2 pipeline with wrapper functions that ease the processing of
  multiple projects

- Some plotting functions for beta diversity, heatmap and differential
  abundance analysis using directly a phyloseq object

- A pipeline for IgASeq analysis

</div>

# Trim, denoise and align your sequences using dada2 pipeline

<div style="text-align: justify">

For convenience this will not be a reproducible example, dada2 takes too
long to compute and knit. This part of the tutorial will present a run
that we performed in house. The rest of the tutorial will be based on
reproducible data.

</div>

## Get the files

<div style="text-align: justify">

We use here a wrapper function that will create a list of three for pair
end :

- forward files

- reverse files

- names of the files

And for single end : -

- a list of files

- names of the files

</div>

``` r

f_list = list_fastq("/home/bigbeast/Documents/tmp/2022-12 CIMMAP run ELISE", pattern = c("R1",
    "R2"), separator = "_", level = 1)

# check that all files are distinct

lapply(f_list, duplicated)

# check that all files exist5
lapply(f_list, file.exists)
random = lapply(f_list, "[", sample(1:255, 30))
# tmp= summarise_fastq(random, cores =30, plot = F )
```

## Check the quality profile

<div style="text-align: justify">

We will use the function `qc_check`. This will take time as the
`plotQualityProfile` isn’t parallelized in dada2. This function will
create two plots (for pair end) of `n` aggregated samples and only one
plot if you are using single end.

</div>

``` r
qc_check(flist, n = 30)
```

## Check the best trimming parameters for your set of fastq

### Test optimal cutting param

``` r
set.seed(1)
tmp = lapply(f_list, "[", sample(1:255, 30))
tmp_list = filt_list(tmp)
filtered = list()

cb = combn(x = (300 - seq(0, 50, by = 10)), m = 2)
# Very long be careful try to test combination of trimmings
for (i in 1:dim(cb)[2]) {
    filtered[[i]] = filter_fastq(tmp, tmp_list, cutting_param = cb[, i], cores = 35,
        trimleft = 35, maxEE = c(3, 4))
}
# extract the data to a df
per = NULL
for (i in 1:15) {
    t = filtered[[i]] %>%
        as.data.frame() %>%
        mutate(per = reads.out/reads.in)
    per = cbind(per, t$per) %>%
        as.data.frame()
    # per$names=rownames(tmp[1])
    print(per)
}

colnames(per) = paste(cb[1, ], cb[2, ])
per$names = rownames(t[1])

# plot

png("transmic and IgA rescue mice percent reads passed depending on cutting parameters.png",
    width = 600, height = 400)
per %>%
    as.data.frame() %>%
    pivot_longer(names_to = "cut", values_to = "per", cols = 1:15) %>%
    ggplot(aes(y = per, x = cut)) + geom_boxplot() + geom_line(aes(group = names,
    col = names)) + labs(y = "percentage passed", x = "Trimming parameters") + theme(axis.text.x = element_text(angle = 90),
    legend.position = "none")


dev.off()
```

### Test optimal maxEE

``` r
cb = combn(x = (7 - seq(0, 5, by = 1)), m = 2)
cb = cb[nrow(cb):1, ]
# Very long be careful try to test combination of trimmings
# registerDoParallel(cl = makeCluster(5))
for (i in 1:dim(cb)[2]) {
    filtered[[i]] = filter_fastq(tmp, tmp_list, cutting_param = c(260, 250), cores = 35,
        trimleft = 35, maxEE = cb[, i])
}
# extract the data to a df
per = NULL
for (i in 1:15) {
    t = filtered[[i]] %>%
        as.data.frame() %>%
        mutate(per = reads.out/reads.in)
    per = cbind(per, t$per) %>%
        as.data.frame()
}

colnames(per) = paste(cb[1, ], cb[2, ])
per$names = rownames(t[1])

# plot

png("transmic and IgA rescue mice percent reads passed depending on errors parameters.png",
    width = 600, height = 400)
per %>%
    as.data.frame() %>%
    pivot_longer(names_to = "cut", values_to = "per", cols = 1:15) %>%
    ggplot(aes(y = per, x = cut)) + geom_boxplot() + geom_line(aes(group = names,
    col = names)) + labs(y = "percentage passed", x = "Trimming parameters") + theme(axis.text.x = element_text(angle = 90),
    legend.position = "none")
dev.off()
```

## Filter and trim the fastq

<div style="text-align: justify">

We now have to remove the bad quality reads and trim the length. You
will find a function to create the list of filtered files and one to
make the filtered files.

</div>

``` r
filt = filt_list(f_list)  # create the list of filtered files

filtered = filterAndTrim(fwd = fwd, filt = filtFs, rev = rv, filt.rev = filtRs, truncLen = c(260,
    240), trimLeft = 25, maxEE = c(3, 5), multithread = 45)
```

# Analyse microbiome system

## Analysis based on phyloseq objects

### Alpha diversity

<div style="text-align: justify">

We will use the **enterotype** data to explore some of the plotting
functions. Let’s start with the beta diversity functions
`beta_diversity` and `beta_dispersion`.

Alpha diversity is an important facet of microbiome analysis. I’ve
created wrapper function to plot either alpha diversity as boxplots or
as line plots. This function is rudimentary but allows to plot quickly
alpha diversity

</div>

``` r
data("enterotype")
alpha_diversity(enterotype, measure = "Shannon", x = "Enterotype", group = "SeqTech",
    plot_type = "boxplot")

alpha_diversity(enterotype, measure = "Shannon", x = "Enterotype", group = "SeqTech",
    plot_type = "line")
```

<img src="ImmuMicrobiome_book_files/figure-gfm/Alpha diversity 1-1.png" style="display: block; margin: auto;" /><img src="ImmuMicrobiome_book_files/figure-gfm/Alpha diversity 1-2.png" style="display: block; margin: auto;" />

<div style="text-align: justify">

These plots are compatible with other aspects of ggplot and tidyverse in
general. You can use facets or stats and even pipes operator to pass
functions before using `alpha_diversity()`. Stats are also implemented
in this function but are very sparse, it can be a good idea to make the
stat on your own wit `stat_compare_means()` or `stat_pvalue_manual()`.

</div>

``` r
alpha_diversity(enterotype, measure = "Shannon", x = "Enterotype", group = "SeqTech",
    plot_type = "boxplot", stat = T)

enterotype %>%
    subset_samples(Enterotype != "NA") %>%
    alpha_diversity(measure = "Shannon", x = "Enterotype", group = "Enterotype",
        plot_type = "boxplot") + stat_compare_means(comparisons = list(c("1", "2"))) +
    facet_grid(~SeqTech)
```

<img src="ImmuMicrobiome_book_files/figure-gfm/Alpha diversity 2-1.png" style="display: block; margin: auto;" /><img src="ImmuMicrobiome_book_files/figure-gfm/Alpha diversity 2-2.png" style="display: block; margin: auto;" />

<div style="text-align: justify">

There is another feature that can become handy at some points, checking
the link between depth and alpha diversity. For this, the function
allows you to specify `check_depth=T`. This will produce two plots:

- one boxplot
- one dotplot with alpha diversity measure as Y and the sum of reads for
  each sample as X

</div>

``` r
data("GlobalPatterns")
GlobalPatterns %>%
    alpha_diversity(x = "SampleType", group = "SampleType", check_depth = T)
```

<img src="ImmuMicrobiome_book_files/figure-gfm/Alpha diversity 3-1.png" style="display: block; margin: auto;" />

### Beta diversity

<div style="text-align: justify">

You will have the choice between `beta_diversity` and `beta_dispersion`
for your beta diversity plotting.

`beta_dispersion` will plot PCoA, NMDS, PCA, DCA, CA and t-SNE for the
moment. This function will plot the two components of your choosing,
**confidence ellipses** and boxplot for each axis and for each group.
Each function will return a plot and a percentage of contribution for
each component.

`beta_diversity` will be removed progressively from this package. This
function is redundant with `beta_dispersion`.

</div>

``` r
beta_diversity(enterotype, dist = "bray", method = "PCoA", group = "SeqTech", permanova = F)

beta_dispersion(enterotype, dist = "bray", method = "PCoA", group = "SeqTech")
```

<img src="ImmuMicrobiome_book_files/figure-gfm/Beta diversity 1-1.png" style="display: block; margin: auto;" /><img src="ImmuMicrobiome_book_files/figure-gfm/Beta diversity 1-2.png" style="display: block; margin: auto;" />

The `beta_dispersion()` function comes with three choices of plotting:

- Boxplot on the side
- Just the multidimensional reduction
- Multidimensional reduction and loadings if the algorithm allows it

NMDS, PCA, CA, DCA are compatible with the `type="arrows"` argument
while PCoA is not. Arrows can be tricky to use with all the taxa in a
phyloseq object. This plotting option is more interesting with sparse
data.

``` r
# Boxplots on the sides
beta_dispersion(enterotype, dist = "bray", method = "PCoA", type = "boxplot", group = "SeqTech",
    stat = "permanova", color_vector = c("#777711", "#117777", "#DD7788"), legend_title = "Sequencing tech",
    lwd = 2, font = 2, draw = "polygon", text = T, y.intersp = 0.7)

# Nothing
beta_dispersion(enterotype, dist = "bray", method = "PCoA", type = "pure", group = "SeqTech",
    stat = "permanova", color_vector = c("#777711", "#117777", "#DD7788"), legend_title = "Sequencing tech",
    lwd = 2, font = 2, draw = "polygon", text = T, y.intersp = 0.7)

# Arrows
beta_dispersion(enterotype, dist = "bray", method = "PCA", type = "arrows", group = "SeqTech",
    stat = "permanova", color_vector = c("#777711", "#117777", "#DD7788"), legend_title = "Sequencing tech",
    lwd = 2, font = 0.1, draw = "polygon", text = T, y.intersp = 0.7)
```

<img src="ImmuMicrobiome_book_files/figure-gfm/Beta diversity 2-1.png" style="display: block; margin: auto;" /><img src="ImmuMicrobiome_book_files/figure-gfm/Beta diversity 2-2.png" style="display: block; margin: auto;" /><img src="ImmuMicrobiome_book_files/figure-gfm/Beta diversity 2-3.png" style="display: block; margin: auto;" />

#### Multiple dimension reductions

<div style="text-align: justify">

There are unconstrained and constrained approaches to beta diversity.
Unconstrained analysis, more common in literature, are either based on:

- dissimilarity matrix (PCoA, NMDS, DCA)
- directly on the raw matrix (PCA, CA)
- t-SNE that can use both dissimilarity matrix or the raw matrix

</div>

##### Unconstrained models

You can play with the parameters of each function like the following
plots :

###### PCoA

``` r
beta_dispersion(enterotype, dist = "bray", method = "PCoA", group = "SeqTech", stat = "permanova",
    color_vector = c("#777711", "#117777", "#DD7788"), legend_title = "Sequencing tech",
    lwd = 2, font = 2, draw = "polygon", text = T, y.intersp = 0.7)
```

<img src="ImmuMicrobiome_book_files/figure-gfm/Beta diversity 3-1.png" style="display: block; margin: auto;" />

###### t-SNE

``` r
beta_dispersion(enterotype, dist = "bray", method = "tsne", group = "SeqTech", color_vector = c("#777711",
    "#117777", "#DD7788"), stat = "permanova", legend_title = "Sequencing tech",
    lwd = 2, font = 2, draw = "polygon", text = T, y.intersp = 0.7, where = "bottom")
```

<img src="ImmuMicrobiome_book_files/figure-gfm/Beta diversity 4-1.png" style="display: block; margin: auto;" />

###### NMDS

``` r
beta_dispersion(enterotype, dist = "bray", method = "NMDS", group = "SeqTech", color_vector = c("#777711",
    "#117777", "#DD7788"), legend_title = "Sequencing tech", lwd = 2, font = 2, draw = "polygon",
    text = T, y.intersp = 0.7, where = "bottom")
#> Run 0 stress 0.1380077 
#> Run 1 stress 0.1536597 
#> Run 2 stress 0.1596971 
#> Run 3 stress 0.1592765 
#> Run 4 stress 0.1499595 
#> Run 5 stress 0.1582167 
#> Run 6 stress 0.1549281 
#> Run 7 stress 0.15729 
#> Run 8 stress 0.1578031 
#> Run 9 stress 0.1521635 
#> Run 10 stress 0.148488 
#> Run 11 stress 0.150579 
#> Run 12 stress 0.1458337 
#> Run 13 stress 0.1539079 
#> Run 14 stress 0.1443015 
#> Run 15 stress 0.1556298 
#> Run 16 stress 0.1508059 
#> Run 17 stress 0.1549227 
#> Run 18 stress 0.1603283 
#> Run 19 stress 0.1614954 
#> Run 20 stress 0.153975 
#> *** Best solution was not repeated -- monoMDS stopping criteria:
#>     18: stress ratio > sratmax
#>      2: scale factor of the gradient < sfgrmin
```

<img src="ImmuMicrobiome_book_files/figure-gfm/Beta diversity 5-1.png" style="display: block; margin: auto;" />

###### DCA

``` r
beta_dispersion(enterotype, dist = "bray", method = "DCA", group = "SeqTech", color_vector = c("#777711",
    "#117777", "#DD7788"), legend_title = "Sequencing tech", lwd = 2, font = 2, draw = "polygon",
    text = T, y.intersp = 0.7, where = "bottom")
```

<img src="ImmuMicrobiome_book_files/figure-gfm/Beta diversity 6-1.png" style="display: block; margin: auto;" />

###### PCA

``` r
beta_dispersion(enterotype, dist = "bray", method = "PCA", group = "SeqTech", color_vector = c("#777711",
    "#117777", "#DD7788"), legend_title = "Sequencing tech without boxplots", lwd = 2,
    font = 2, draw = "polygon", text = T, y.intersp = 0.7, where = "bottom")
```

<img src="ImmuMicrobiome_book_files/figure-gfm/Beta diversity 7-1.png" style="display: block; margin: auto;" />

###### CA

``` r
beta_dispersion(enterotype, dist = "bray", method = "CA", group = "SeqTech", color_vector = c("#777711",
    "#117777", "#DD7788"), legend_title = "Sequencing tech without boxplots", lwd = 2,
    font = 2, draw = "polygon", text = T, y.intersp = 0.7, where = "bottom")
```

<img src="ImmuMicrobiome_book_files/figure-gfm/Beta diversity 8-1.png" style="display: block; margin: auto;" />
\##### {-}

#### Constrained models

Constrained analysis are multiple regression of the unconstrained
approaches. They allow to fit the ecological data or sample data to the
community data. Alternatively, you can fit other omic data to your beta
diversity analysis. We have two different classes of constrained
analysis :  

- `rda` is a regression of a PCA
- `cca` is a regression of a CA
- `dbRDA` is a RDA based on a dissimilarity matrix  
  The function will return a plot, same as `beta_dispersion`, and also
  the result of the constrained analysis.

Vegan authors have implemented a very nice feature if you want to make
your own representation instead of using mine. You can set use
`scores(res, tidy=T)` to return a dataframe compatible with ggplot2.

``` r
mod = "SeqTech+Gender+Nationality+Age+ClinicalStatus"
```

###### CCA

``` r
res = constrained_beta_dispersion(enterotype, model = mod, group = "Enterotype",
    method = "CCA", text = T, color_vector = tol21rainbow)
```

<img src="ImmuMicrobiome_book_files/figure-gfm/Beta diversity 9-1.png" style="display: block; margin: auto;" />

###### RDA

``` r
res = constrained_beta_dispersion(enterotype, model = mod, group = "Enterotype",
    method = "RDA", text = T, color_vector = tol21rainbow)
```

<img src="ImmuMicrobiome_book_files/figure-gfm/Beta diversity 10-1.png" style="display: block; margin: auto;" />

###### dbRDA

``` r
res = constrained_beta_dispersion(enterotype, model = mod, group = "Enterotype",
    method = "dbRDA", text = T, color_vector = tol21rainbow)
```

<img src="ImmuMicrobiome_book_files/figure-gfm/Beta diversity 11-1.png" style="display: block; margin: auto;" />

You can then use `anova` statistic test to decipher which factor or
feature is significantly associated with your beta diversity.

``` r
anova(res, by = "terms")
#> Permutation test for capscale under reduced model
#> Terms added sequentially (first to last)
#> Permutation: free
#> Number of permutations: 999
#> 
#> Model: capscale(formula = as(otu_table(reverseASV(physeq)), "matrix") ~ SeqTech + Gender + Nationality + Age + ClinicalStatus, data = df, distance = dist, na.action = na.exclude)
#>                Df SumOfSqs      F Pr(>F)    
#> Gender          1  0.06650 0.7390  0.639    
#> Nationality     5  1.06616 2.3695  0.005 ** 
#> Age             1  0.43786 4.8657  0.001 ***
#> ClinicalStatus  3  0.30763 1.1395  0.291    
#> Residual       26  2.33975                  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Heatmap based on the phyloseq object

<div style="text-align: justify">

This function will perform a top taxa at a rank of your choosing and
create a heatmap with annotations. For now only one annotation is
supported.  
The clustering is made using the `hclust`function and a Ward.D2 method.
You can define the distance matrix that you want to use. It is important
to note that the distance is made before any trimming of the data. This
means that the distance matrix is made at the ASV/OTU level before doing
the rank merging and topping, so the clustering will represent your
“true” data instead of a modified dataset. In my personnal opinion, this
is the only way of performing clustering or any type of generalized
analysis : clustering, reducing or else on the ASV/OTU/MAG/… level and
then aggregates for plotting.

If you really want to cluster on other taxonomic level you should use
`tax_glom` before using this function.

``` r
phylo_heatmap(enterotype, top = 30, labels = "SeqTech", taxa_rank = "Genus", factor_to_plot = "Enterotype",
    split = 3, distance = "bray")
#> [1] "No phylogenetic tree in this phyloseq object, bray-curtis distance selected."
#> [1] " not reversed"
```

<img src="ImmuMicrobiome_book_files/figure-gfm/Phylo heatmap 1-1.png" style="display: block; margin: auto;" />

</div>

## Usage of various functions

<div style="text-align: justify">

The idea behind these functions is : creating a more automatic pipeline
enabling filtering and subseting of the `phyloseq` object without having
to perform a temporary `phyloseq` object and a temporary distance
object. With these function you can directly use the “pipe” introduced
by `magrittr`.  
So you can use a `subset_samples` like the following and automatically
plot your beta diversity without adding too much code.

``` r
enterotype %>%
    subset_samples(SeqTech == "Illumina") %>%
    beta_dispersion(group = "Enterotype", color_vector = c("#777711", "#117777",
        "#DD7788"), legend_title = "Enterotypes \n with bray", lwd = 2, font = 2,
        draw = "lines", text = T, stat = "permanova", y.intersp = 0.7, where = "bottomleft",
        cex = 3)
```

<img src="ImmuMicrobiome_book_files/figure-gfm/Phylo heatmap 2-1.png" style="display: block; margin: auto;" />

Additionally, if you want to go further you can also serialized the code
with a `for` loop or a `lapply` or even a parallel approach using
`mclapply` or `doParallel`.

``` r
layout(matrix(c(1, 2, 3), nrow = 1, ncol = 3, byrow = TRUE))
for (i in c("Illumina", "Sanger", "Pyro454")) {
    enterotype %>%
        subset_samples(SeqTech == i & !is.na(Enterotype)) %>%
        beta_dispersion(group = "Enterotype", color_vector = c("#777711", "#117777",
            "#DD7788"), legend_title = "Enterotypes \n with bray", type = "pure",
            lwd = 2, font = 2, draw = "lines", text = T, stat = "none", y.intersp = 0.7,
            where = "bottomleft", cex = 3)
}
```

<img src="ImmuMicrobiome_book_files/figure-gfm/unnamed-chunk-8-1.png" style="display: block; margin: auto;" /><img src="ImmuMicrobiome_book_files/figure-gfm/unnamed-chunk-8-2.png" style="display: block; margin: auto;" /><img src="ImmuMicrobiome_book_files/figure-gfm/unnamed-chunk-8-3.png" style="display: block; margin: auto;" />

</div>

## Differential abundance testing

Differential abundance testing is quite complicated because of the
number of different statistical approaches existing in the literature. I
personally use ALDEx2 and SIAMCAT a lot. The first because it is usually
highly regarded and the second for the overall easiness of the package.
For now `differential_abundance` implement only `ALDEx2` approach.

The function will return two datasets (all taxa and one with only the
significant taxa) and two plots representing the taxa significantly
represented in one of the two conditions. This function only covers two
by two analysis.

``` r
data("GlobalPatterns")
tmp = GlobalPatterns %>%
    subset_samples(SampleType == "Feces" | SampleType == "Soil") %>%
    subset_taxa(Genus != "NA") %>%
    tax_glom("Genus")
taxa_names(tmp) = paste0(tax_table(tmp)[, "Genus"], 1:length(tax_table(tmp)[, "Genus"]))

res = tmp %>%
    differential_abundance(group = "SampleType", col1 = "brown", col2 = "darkgreen",
        plot = T)
```

### Results

``` r
head(res$all_features)
#>                                 we.ep     we.eBH      wi.ep    wi.eBH
#> Nitrosopumilus3           0.038020857 0.12553457 0.07321429 0.1739762
#> CandidatusNitrososphaera4 0.080883984 0.21608822 0.05714286 0.1534536
#> Methanocorpusculum15      0.393517314 0.52330393 0.50937500 0.6131738
#> Methanobacterium17        0.192887050 0.32809281 0.22232143 0.3464315
#> Methanobrevibacter19      0.005763966 0.04577142 0.05714286 0.1534536
#> Propionibacterium20       0.575347882 0.69177092 0.69397321 0.7805085
#>                               rab.all rab.win.Feces rab.win.Soil    diff.btw
#> Nitrosopumilus3           -0.03578414     2.1640049  -3.54657967  -6.0725458
#> CandidatusNitrososphaera4  2.94274075     1.7275926   6.88742882   5.5007228
#> Methanocorpusculum15      -0.74114455    -1.6164717  -0.03037559   1.3795403
#> Methanobacterium17        -1.58879597     0.6563627  -3.78683212  -4.4195807
#> Methanobrevibacter19       5.47119787     9.1212605  -3.52760211 -13.2996558
#> Propionibacterium20       -2.23757862    -1.6605794  -3.39678481  -0.9666516
#>                           diff.win     effect      overlap
#> Nitrosopumilus3           3.082815 -1.8944361 0.0155616731
#> CandidatusNitrososphaera4 3.935999  1.8209332 0.0003653998
#> Methanocorpusculum15      3.092592  0.4041233 0.2953373376
#> Methanobacterium17        4.154777 -0.8932739 0.1398980468
#> Methanobrevibacter19      4.791323 -2.9130433 0.0003653998
#> Propionibacterium20       5.166144 -0.1619715 0.4322918226
```

### Barplot results

``` r
res$barplot
```

<img src="ImmuMicrobiome_book_files/figure-gfm/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

### Volcano plot results

``` r
res$volcano
```

<img src="ImmuMicrobiome_book_files/figure-gfm/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

## Non phyloseq object oriented functions

::: {style=“text-align: justify”} To enforce same plotting capacity on
non phyloseq objects, the `plot_reduction` and
`plot_constrained_reduction` have been created. These functions needs a
data.frame containing categorical data and the numeric values to make
the reduction. These functions are perfect to analyse metabolomic data,
flow cytometry data or else. Usage remains the same, with the same
arguments than in `beta_dispersion` and `constrained_beta_dispersion`.

### Reduction of dimensions on non phyloseq objects

``` r
plot_reduction()

plot_constrained_reduction()
```

The first function is an equivalent of and the second one is an
equivalent of . Make sure to use data.frame with samples as rows.

``` r
data(metabolomic)
plot_reduction(mat = metabolomic, clinical_data = 1:4, group = "birth_type", method = "CA",
    type = "pure", stat = "permanova")

plot_reduction(mat = metabolomic, clinical_data = 1:4, group = "birth_type", method = "PCA",
    type = "pure", stat = "permanova")
plot_reduction(mat = metabolomic, clinical_data = 1:4, group = "birth_type", method = "PCoA",
    dist = "euclidean", type = "pure", stat = "permanova")
```

<img src="ImmuMicrobiome_book_files/figure-gfm/reduction2-1.png" style="display: block; margin: auto;" /><img src="ImmuMicrobiome_book_files/figure-gfm/reduction2-2.png" style="display: block; margin: auto;" /><img src="ImmuMicrobiome_book_files/figure-gfm/reduction2-3.png" style="display: block; margin: auto;" />

### PLS-DA

<div style="text-align: justify">

Partial Least Square Discriminant Analysis are widely used, especially
in metabolomics. This particular statistical approach is not yet
supported in `constrained_beta_dispersion`and
`plot_constrained_reduction` but will be in the late versions. For now,
this function will use the results from the function and use ggplot2 to
make the graph. You can customize the graph using classic `ggplot2`
functions.

``` r
results = mixOmics::plsda(X = metabolomic[, 5:217], Y = metabolomic$birth_type)

# basic graphic

plot_plsda(results, top.loads = 20, color_vector = c("brown", "orange"))

# with a bit of improvement

plot_plsda(results, top.loads = 20, color_vector = c("brown", "orange")) + guides(fill = guide_legend("Birth route")) +
    theme_bw() + theme(legend.position = "bottom", axis.title = element_text(size = 15,
    face = "bold"), panel.grid = element_blank(), legend.text = element_text(size = 15),
    legend.title = element_text(size = 20, face = "bold"), axis.text = element_text(size = 15))
```

<img src="ImmuMicrobiome_book_files/figure-gfm/plsda1-1.png" style="display: block; margin: auto;" /><img src="ImmuMicrobiome_book_files/figure-gfm/plsda1-2.png" style="display: block; margin: auto;" />

</div>

### Analysis of variance

<div style="text-align: justify">

Here is an important part of any analysis : finding which factors
explains the most variance. The principle is that we will get the
variance for each features and the variance for each group. The group
variance will be summed for the given factor and then we will use
1-variance.factor/variance.total, which will give us the variance for
each factors.

On top of the variance, the function will calculate p.values using
kruskal tests or wilcoxon test depending of the number of factors and a
FDR correction. For now the function only supports non parametric tests.

The function can take quite some time if you have a lot of samples,
features and factors, you can use multi-threading here.

``` r
var = metabolomic %>%
    dplyr::select(!child_id) %>%
    calculate_variance(clinical_data = 1:3, cores = 1)

lapply(var, head, 5)
#> $variance
#> # A tibble: 5 × 5
#>   features               variable   birth_type breastfeeding          sex
#>   <chr>                  <fct>           <dbl>         <dbl>        <dbl>
#> 1 1,5-Naphthalenediamine var.tot  0.0000000137  0.0000000137 0.0000000137
#> 2 1,7-Dimethyluric acid  var.tot  0.0000000855  0.0000000855 0.0000000855
#> 3 1-Methyladenine        var.tot  0.0000000440  0.0000000440 0.0000000440
#> 4 1-Methylguanine        var.tot  0.000000345   0.000000345  0.000000345 
#> 5 1-Vinylimidazole       var.tot  0.000000394   0.000000394  0.000000394 
#> 
#> $p.value
#> # A tibble: 5 × 4
#>   features                    p p.adj factor    
#>   <chr>                   <dbl> <dbl> <chr>     
#> 1 1,5-Naphthalenediamine 0.0993 0.414 birth_type
#> 2 1,7-Dimethyluric acid  0.609  0.786 birth_type
#> 3 1-Methyladenine        0.0555 0.404 birth_type
#> 4 1-Methylguanine        0.266  0.596 birth_type
#> 5 1-Vinylimidazole       0.64   0.807 birth_type
#> 
#> $mean.feat
#>                 features    mean.feat
#> 1 1,5-Naphthalenediamine 0.0001591865
#> 2  1,7-Dimethyluric acid 0.0002805524
#> 3        1-Methyladenine 0.0001991781
#> 4        1-Methylguanine 0.0003367005
#> 5       1-Vinylimidazole 0.0001454276
#> 
#> $all
#> $all$birth_type
#> # A tibble: 213 × 5
#>    features                        var.tot mean.feat  var.grp var.exp
#>    <chr>                             <dbl>     <dbl>    <dbl>   <dbl>
#>  1 1,5-Naphthalenediamine         1.37e- 8 0.000159  1.27e- 8  0.0768
#>  2 1,7-Dimethyluric acid          8.55e- 8 0.000281  8.44e- 8  0.0130
#>  3 1-Methyladenine                4.40e- 8 0.000199  4.27e- 8  0.0313
#>  4 1-Methylguanine                3.45e- 7 0.000337  3.37e- 7  0.0238
#>  5 1-Vinylimidazole               3.94e- 7 0.000145  3.88e- 7  0.0158
#>  6 17α-Hydroxyprogesterone        6.76e- 8 0.000156  6.64e- 8  0.0169
#>  7 1H-indene-3-carboxamide        8.47e-10 0.0000441 8.12e-10  0.0411
#>  8 2,3-pyridinecarboxylic acid    1.25e- 7 0.000426  1.23e- 7  0.0125
#>  9 2-(hydroxymethyl)butanoic acid 2.06e- 6 0.000646  2.02e- 6  0.0150
#> 10 2-Fucosyllactose               1.50e- 5 0.00142   1.43e- 5  0.0482
#> # ℹ 203 more rows
#> 
#> $all$breastfeeding
#> # A tibble: 213 × 5
#>    features                        var.tot mean.feat  var.grp var.exp
#>    <chr>                             <dbl>     <dbl>    <dbl>   <dbl>
#>  1 1,5-Naphthalenediamine         1.37e- 8 0.000159  1.31e- 8  0.0464
#>  2 1,7-Dimethyluric acid          8.55e- 8 0.000281  8.44e- 8  0.0134
#>  3 1-Methyladenine                4.40e- 8 0.000199  4.35e- 8  0.0126
#>  4 1-Methylguanine                3.45e- 7 0.000337  3.40e- 7  0.0154
#>  5 1-Vinylimidazole               3.94e- 7 0.000145  3.89e- 7  0.0117
#>  6 17α-Hydroxyprogesterone        6.76e- 8 0.000156  6.56e- 8  0.0291
#>  7 1H-indene-3-carboxamide        8.47e-10 0.0000441 8.13e-10  0.0404
#>  8 2,3-pyridinecarboxylic acid    1.25e- 7 0.000426  1.23e- 7  0.0108
#>  9 2-(hydroxymethyl)butanoic acid 2.06e- 6 0.000646  2.03e- 6  0.0113
#> 10 2-Fucosyllactose               1.50e- 5 0.00142   1.48e- 5  0.0139
#> # ℹ 203 more rows
#> 
#> $all$sex
#> # A tibble: 213 × 5
#>    features                        var.tot mean.feat  var.grp var.exp
#>    <chr>                             <dbl>     <dbl>    <dbl>   <dbl>
#>  1 1,5-Naphthalenediamine         1.37e- 8 0.000159  1.33e- 8  0.0303
#>  2 1,7-Dimethyluric acid          8.55e- 8 0.000281  8.46e- 8  0.0106
#>  3 1-Methyladenine                4.40e- 8 0.000199  3.92e- 8  0.109 
#>  4 1-Methylguanine                3.45e- 7 0.000337  3.39e- 7  0.0192
#>  5 1-Vinylimidazole               3.94e- 7 0.000145  3.83e- 7  0.0270
#>  6 17α-Hydroxyprogesterone        6.76e- 8 0.000156  6.64e- 8  0.0174
#>  7 1H-indene-3-carboxamide        8.47e-10 0.0000441 8.27e-10  0.0230
#>  8 2,3-pyridinecarboxylic acid    1.25e- 7 0.000426  1.23e- 7  0.0102
#>  9 2-(hydroxymethyl)butanoic acid 2.06e- 6 0.000646  2.01e- 6  0.0224
#> 10 2-Fucosyllactose               1.50e- 5 0.00142   1.48e- 5  0.0172
#> # ℹ 203 more rows
```

The function will produce a list of data.frame, one for each factor in
your dataset. In each data.frame you will find:

- variance: a data.frame with total variance, the group variance, the
  explained variance and the mean quantity of a given features
  - “var.tot” aka total variance for the factor, all NA value for the
    given factor are removed, so the total variance can differ between
    factors
  - “var.grp” the variance explained by the factors
  - “var.exp” explained variance as described by
    $$explained variance = 1 - group variance / total variance$$
  - “mean.feat” the mean value for a given feature
- p.value : a data.frame with the p values and p values corrected with
  FDR for a given feature and a given factor. Currently, the function
  support only non parametric tests.
- mean.feat: the mean value for a given feature

</div>

#### Plot variance

<div style="text-align: justify">

Now that the variance is calculated for all your factors, you need to
visualise it. I’ve come up with different ways to visualise the variance
of a given dataset:

- `plot_all_variance`
- `plot_xy_variance`
- `circular_variance_plot`

</div>

##### Plot all variance

<div style="text-align: justify">

`plot_all_variance` will give a boxplot graph with a dot for each
features, the size is determined by the mean.feat value for a given
feature and the color is based on the p.values.

``` r
plot_all_variance(var, col = c("brown", "darkgreen", "grey"))
```

<img src="ImmuMicrobiome_book_files/figure-gfm/variance 2-1.png" style="display: block; margin: auto;" />
Alternatively, the function accepts the `heatmap`argument, this will
produce a basic heatmap. However, this part of the function is still
under development at the moment. In the end the function will provide
the possibility to clusterise features.

``` r
plot_all_variance(var, plot_type = "heatmap")
```

<img src="ImmuMicrobiome_book_files/figure-gfm/variance 3-1.png" style="display: block; margin: auto;" />

</div>

##### Plot variance two by two

<div style="text-align: justify">

`plot_xy_variance` will return a dot plot with factor 1 as X and factor
2 as Y. The p-values are displayed as following: “Both factors”, “factor
1”, “factor 2”, “None”.

``` r
plot_xy_variance(var, x = "birth_type", y = "breastfeeding", corrected = F)
```

<img src="ImmuMicrobiome_book_files/figure-gfm/variance 4-1.png" style="display: block; margin: auto;" />

</div>

##### Plot variance as circular plot

Another way to visualise the variance in the dataset is to use the
`circular_variance_plot` function. This function is not really suitable
when many features are present.

``` r

circular_variance_plot(var)

metabolomic %>%
    select(2:20) %>%
    calculate_variance(clinical_data = 1:3, cores = 1) %>%
    circular_variance_plot(adjust = 0.3) + ylim(c(-0.05, 0.5))
```

<img src="ImmuMicrobiome_book_files/figure-gfm/variance 5-1.png" style="display: block; margin: auto;" /><img src="ImmuMicrobiome_book_files/figure-gfm/variance 5-2.png" style="display: block; margin: auto;" />

### Correlations

<div style="text-align: justify">

One last aspect of classical biostatisics is correlations. It can be
extremely handy to verify correlation between taxa and metabolomics for
example. Correlation can be presented as dotplot with a regression line
or as heatmaps/correlograms. Good packages for plotting correlogramms
already exists on R, I’ve decided to implement a heatmap correlation
graph.  
For thsi example we will take only 30 random features to ensure plotting
capacity.

</div>

``` r
# take randomly 30 features from the metabolomic dataset
rng = sample(x = 5:217, 30)
rng2 = sample(x = 5:217, 30)
plot_corr_heatmap(X = metabolomic[, rng], Y = metabolomic[, rng2], ratio = 0.5) +
    theme(text = element_text(size = 10))
```

<img src="ImmuMicrobiome_book_files/figure-gfm/correlation 1-1.png" style="display: block; margin: auto;" />

<div style="text-align: justify">

By default the function will keep every features, you can specify a
cutoff of correlation factor, here we will keep only the features were
at least one as a correlation factor \> 0.7 and with a p.adj \< 0.05.
The function will only return features checking those two criteria, for
now you can’t chose between corrected and non corrected p values.  
You can specify two different matrices, for example the taxa table and
the features to find correlations between these, just make sure you’re
data are in the correct order because the function assumes that the
samples are matching between the two matrices.

</div>

``` r
plot_corr_heatmap(X = metabolomic[, rng], Y = metabolomic[, rng2], cutoff = 0.7)
```

<img src="ImmuMicrobiome_book_files/figure-gfm/correlation 2-1.png" style="display: block; margin: auto;" />

Finally, you can also chose to cluster the features both for x and y
axis.

``` r
plot_corr_heatmap(X = metabolomic[, rng], Y = metabolomic[, rng2], cutoff = 0.4,
    cluster = T, ratio = 0.5)
```

<img src="ImmuMicrobiome_book_files/figure-gfm/correlation 3-1.png" style="display: block; margin: auto;" />

### Volcano and other plots to represent enrichment

<div style="text-align: justify">

Finally, we finish with enrichment/differential abundance analysis for
non phyloseq objects. This is a big part of biostat analysis, find the
features/taxa differently abundant between two groups of patient/mice.
There are a lot of discussion about how to normalize data in the
microbiomics, this might be the same case for metabolomics and
metagenomics. Here the function use a classic Mann-Whitney or Student
test with a FDR correction. So this is the most basic approach. You can
chose to normalize the data before giving it to the function, this will
be on you to chose which approach is better.

#### Volcano plots

We already have a function producing volcano plots based on the ALDeX2
package but not for non phyloseq objects.

``` r
fold = fold_change(clinical_data = 1:4, metabolomic, cores = 1)
plot_volcano(fold$birth_type)
#> Warning: ggrepel: 6 unlabeled data points (too many overlaps). Consider
#> increasing max.overlaps
```

<img src="ImmuMicrobiome_book_files/figure-gfm/fold change-1.png" style="display: block; margin: auto;" />

#### Ternary plots

The function is not yet created but this following code will create a
ternary plot. Bare in mind that this function will accept only three
groups, no more no less. For two groups use the volcano function, for
more than three groups well you nothing exists yet.

</div>

# IgA seq analysis pipeline

<div style="text-align: justify">

You’ve performed all your sortings and sequencing, you now have three
samples coming from a single individual.  
We created a pipeline analysis where you will use the **neg1** (9/10 of
the neg fraction) and **neg2** (1/10 of the neg fraction) dispersion
(centered and reduced) to create a normal dispersion, using a Z approach
we will have a Z score for the based on the standard deviation of .

Here we can see what the analysis will look like:

<div class="figure" style="text-align: center">

<img src="../Documentation_ImmuMicrobiome/Z test.png" alt="Scheme of Z test approach" width="2700" />
<p class="caption">
Scheme of Z test approach
</p>

</div>

- The technical dispersion is assessed using `neg1/neg2`for the log2
  ratio and the `neg1*neg2`abundance for log10 (black dot)

- The biological dispersion is assessed `pos/neg1` for the log2 ratio
  and the `pos*neg1` abundance for log10 (orange dot)

We will then take windows of X ASV (for example 20) to create n Gaussian
curve and n …  
The pipeline is based on three main functions that will call for other
functions : seq_table, slide_z and collapse_IgAseq.

</div>

## Make the seq_table list

<div style="text-align: justify">

First the seq_table function will take your phyloseq object and
transform it to a list of data frames, one data frame for each samples
coming from a single individual. For this function you will need to give
physeq, sample_name corresponding to the name identifying the individual
from which the sorted samples came, sorting_names the column where we
find the samples such as : “sample1_pos”, “sample1_neg1”,
“sample1_neg2”. Then the cols_to_keep that need to stay for now on
“all”, it will collapse your sample_data in one column to allow you get
it back later on.

The function will tell you if there is some samples are alone, if you
have duplicated samples you need to sort them out or the rest of the
pipeline will block. You can take a look at the architecture of the new
object, you will find the ASV sequence as rownames, the taxonomy
collapsed with “\#” separator, the three samples having their own column
and the sample_data collapse using also “\#” separator. The function
will only take the ASV that are present in the samples.

``` r
data("igaseq")
igaseq = transform_sample_counts(igaseq, function(x) x/sum(x))
# sample_names(igaseq)= sample_data(igaseq)$sort_type

seq.tab = seq_table(igaseq, sample_name = "sample_origin", sorting_names = "sample_sort",
    cols_to_keep = "all")
#> sample MO308 is alone 
#> One sample belonging to C1636 has no reads
#> sample C1232 is alone 
#> sample C1283 is alone 
#> sample T2648 is alone
```

``` r

knitr::kable(tail(seq.tab$C1101), rownames = F, booktabs = T)
```

| FALSE                                                                                                                                                                                                                                                                                                                                                                                                                                               | taxonomy                                                                                                     | C1101_pos | C1101_neg2 | C1101_neg1 | sample_id | new                                                |
|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------|----------:|-----------:|-----------:|:----------|:---------------------------------------------------|
| AGTGGGGAATATTGGGCAATGGGGGAAACCCTGACCCAGCAACGCCGCGTGAAGGAAGAAGGCCTTCGGGTTGTAAACTTCTTTTACCAGGGACGAAGGACGTGACGGTACCTGGAGAAAAAGCAACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTTGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGAGATGCAAGTTGGGAGTGAAATCCATGGGCTCAACCCATGAACTGCTCTCAAAACTGTATCCCTTGAGTATCGGAGAGGCAAGCGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGGCTTGCTGGACGACAACTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACAGGATTA                       | Bacteria#Firmicutes#Clostridia#Oscillospirales#Oscillospiraceae#UCG-005#NA#UCG-005_372                       | 0.0000000 |  0.0011692 |  0.0008873 | C1101     | C1101_pos#C1101#FAFR101#M#breastfed#Vaginal#infant |
| AGTGGGGGATATTGCACAATGGGGGAAACCCTGATGCAGCAACGCCGCGTGAGGGAAGAAGGTTTTCGGATTGTAAACCTCTGTTCTTAGTGACGATAATGACGGTAGCTAAGGAGAAAGCTCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGAGCGAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGTAGGCGGCGAGGCAAGTCAGGCGTGAAATCTATGGGCTTAACCCATAAACTGCGCTTGAAACTGTCTTGCTTGAGTGAAGTAGAGGTAGGCGGAATTCCCGGTGTAGCGGTGAAATGCGTAGAGATCGGGAGGAACACCAGTGGCGAAGGCGGCCTACTGGGCTTTAACTGACGCTGAAGCACGAAAGCATGGGTAGCAAACAGGATTA                          | Bacteria#Firmicutes#Clostridia#Oscillospirales#Ruminococcaceae#Incertae Sedis#NA#Incertae Sedis_437          | 0.0000000 |  0.0017987 |  0.0023662 | C1101     | C1101_pos#C1101#FAFR101#M#breastfed#Vaginal#infant |
| AGTCGGGAATATTGCGCAATGGAGGAAACTCTGACGCAGTGACGCCGCGTATGGGAAGAAGGTTTTCGGATTGTAAACCATTTTAGACAAGGAAGAAACAAGACAGTACTTGTAGAATAAGCTCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGAGCAAGCGTTATCCGGATTTATTGGGTGTAAAGGGTGCGTAGACGGGAAGGTAAGTTAGTTGTGAAATCCCTCGGCTCAACTGAGGAACTGCGACTAAAACTGCTTTTCTTGAGTGCTGGAGAGGAAAGTGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGACAGTAACTGACGTTGAGGCACGAAAGTGTGGGGAGCAAACAGGATTA                         | Bacteria#Firmicutes#Clostridia#Clostridia UCG-014#NA#NA#NA#NA_486                                            | 0.0012804 |  0.0009893 |  0.0000000 | C1101     | C1101_pos#C1101#FAFR101#M#breastfed#Vaginal#infant |
| AGTGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCAACGCCGCGTGAGTGAAGAAGTATTTCGGTATGTAAAGCTCTATCAGCAGGGAAGAAAATGACGGTACCTGACTAAGAAGCACCGGCTAAATACGTGCCAGCAGCCGCGGTAATACGTATGGTGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGCAGGCGGTCTGGCAAGTCTGATGTGAAAGCCCGGGGCTCAACCCCGGGACTGCATTGGAAACTGTCAGACTAGAGTGTCGGAGAGGTAAGTGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGGCTTACTGGACGATAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAGGATTA                          | Bacteria#Firmicutes#Clostridia#Lachnospirales#Lachnospiraceae#NA#NA#NA_635                                   | 0.0000000 |  0.0045867 |  0.0040422 | C1101     | C1101_pos#C1101#FAFR101#M#breastfed#Vaginal#infant |
| AGTGGGGAATATTGCACAATGGAGGAAACTCTGATGCAGCGATGCCGCGTGAGGGAAGAAGGTTTTAGGATTGTAAACCTCTGTCTTCAGGGACGAAAAAAAAGACGGTACCTGAGGAGGAAGCTCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGAGCGAGCGTTGTCCGGAATTACTGGGTGTAAAGGGAGCGTAGGCGGGATCGCAAGTCAGATGTGAAAACTATGGGCTTAACCCATAAACTGCATTTGAAACTGTGGTTCTTGAGTGAAGTAGAGGTAAGCGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACATCAGTGGCGAAGGCGGCTTACTGGGCTTTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAGGATTA                       | Bacteria#Firmicutes#Clostridia#Oscillospirales#Ruminococcaceae#Ruminococcus#bicirculans#Ruminococcus_785     | 0.0021668 |  0.0041371 |  0.0020704 | C1101     | C1101_pos#C1101#FAFR101#M#breastfed#Vaginal#infant |
| AGTGGGGAATCTTCCGCAATGGGCGAAAGCCTGACGGAGCAACGCCGCGTGAGTGATGACGGCCTTCGGGTTGTAAAGCTCTGTGATCGGGGACGAACGGTCTGTAAGCTAATATCTTATGGAAGTGACGGTACCCGAATAGCAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGCGGCTTTCTAAGTCCATCTTAAAAGTGCGGGGCTTAACCCCGTGATGGGATGGAAACTGGAAAGCTGGAGTATCGGAGAGGAAAGTGGAATTCCTAGTGTAGCGGTGAAATGCGTAGAGATTAGGAAGAACACCGGTGGCGAAGGCGACTTTCTGGACGACAACTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACAGGATTA | Bacteria#Firmicutes#Negativicutes#Veillonellales-Selenomonadales#Veillonellaceae#Dialister#NA#Dialister_1138 | 0.0000000 |  0.0046767 |  0.0054225 | C1101     | C1101_pos#C1101#FAFR101#M#breastfed#Vaginal#infant |

The colnames will have the sample_names as given in the otu_table(),
check that they correspond to your pos, neg1 and neg2.

Run the main function Now we will run the main function : slide_z. This
function will take you seq_table object and run the Z function for each
samples. If you are running this function for the first time use plot=T,
if you already made the plots let it as FALSE it will make the loop a
lot faster.

</div>

## Deal with the zero values

<div style="text-align: justify">

In this approach we will use log2 ratio and log10 of abundances. As you
know log2(0/x) or log2(x/0) can’t be performed, so we need a way to deal
with the zeros. We came up with two approach :

- remove all ASV with a zero value

- replace 0 by a random number between 0 and the min value found in one
  of the three samples

- replace 0 by the minimum count in each samples, this probably the
  worst thing to do

In sample with few ASV, i.e meconiums, I strongly recommend using the
random_generation while in adults samples the no_zero approach is
performing better. For more complex samples the decision is up to you,
you can use the function neg_dipsersion to visualize the two outcomes.

``` r
neg_dispersion(seq.tab[[3]], positive_sorted_sample = "pos", negative_sorted_sample = "neg1",
    second_negative_sample = "neg2", type = "superposed")

neg_dispersion(seq.tab[[3]], positive_sorted_sample = "pos", negative_sorted_sample = "neg1",
    second_negative_sample = "neg2", type = "facet all three")
```

<div class="figure" style="text-align: center">

<img src="ImmuMicrobiome_book_files/figure-gfm/unnamed-chunk-16-1.png" alt="Verify the dispersion of the negative fractions"  /><img src="ImmuMicrobiome_book_files/figure-gfm/unnamed-chunk-16-2.png" alt="Verify the dispersion of the negative fractions"  />
<p class="caption">
Verify the dispersion of the negative fractions
</p>

</div>

``` r
# run the following if you want every samples, make sure to change the number
# of cores pdf('test.pdf', width = 10, height = 7) mclapply(seq,
# neg_dispersion, mc.cores = 6, type='superposed') dev.off()
```

</div>

## Run the slide_z function

<div style="text-align: justify">

The first step will be to create log2 ratios and log10 abundance (log10
abundance of **pos** \* **neg1** for example) for each ASV between the
**pos** and the **neg1** and between the ***neg1*** and **neg2** . This
will be done by the function `log_ratio` called by the `slide_z`
function.  
The function will also create an ellipse of confidence interval of your
choosing, default being `confidence_interval= c(0.95,0.99, 0.999)`, this
is done by the function `ellipse_me` also called by the `slide_z`
function.

The output will be a list of S4 objects containing the following slots :

- *ig_seq_all* : containing all the samples and all the ASV

- *ig_up* : containing all the samples and all the ASV that
  significantly enriched in the **IgA positive fraction**

- *ig_down* : containing all the samples and all the ASV that are
  significantly enriched in the **IgA negative fraction**

In each slot you will find the following columns :

- *taxonomy*

- *sample_id*

- *new* : the sample_data collapsed

- *pos* : positive fraction abundance

- *neg1* : negative (9/10) fraction abundance

- *neg2* : negative (1/10) fraction abundance

- *log10_abundance* = log10(pos \* neg1)

- *log2_ratio* = log2(pos / neg1)

- *log10_neg_abundance* = log10(neg1 \* neg2)

- *log2_neg_ratio* = log2(neg1 / neg2)

- *taxa* = tax_table collapsed

- *SlideNorm* = the normalized dispersion for each ASV

- *score* = IgAseq score

- *ellipse_level* = the level of confidence

``` r
IgA_seq = list()

for (i in names(seq.tab)) {
    # print(i)
    IgA_seq[[i]] = slide_z(seq.tab[[i]], positive_sorted_sample = "pos", negative_sorted_sample = "neg1",
        second_negative_sample = "neg2", deltaX = 30, slide_version = "slide_z_modern",
        alpha = 0.05, plot = F, zero_treatment = "random generation")
}
```

Just for the demo we will make the plots with `plot= T`. This will plot
.png and .hmtl `plotly`files. This will make the loop slower but it will
allow you to verify the distribution of your samples and potentially
enable you to remove outliers.

``` r
tmp = slide_z(seq.tab[[1]], positive_sorted_sample = "pos", negative_sorted_sample = "neg1",
    second_negative_sample = "neg2", deltaX = 30, slide_version = "slide_z_modern",
    alpha = 0.05, plot = T, zero_treatment = "random generation")
```

<div class="figure" style="text-align: center">

<img src="Igaseq plots/png/MO101_random generation_.png" alt="Example of the plots output, a .png and a .html are saved" width="2100" />
<p class="caption">
Example of the plots output, a .png and a .html are saved
</p>

</div>

</div>

## Collapse the list of IgAseq into a single S4 object

<div style="text-align: justify">

We now need to collapse the list of IgA_seq objects into a single one,
just use `collapse_IgAseq` and separate the columns that were pasted
together.

``` r
IgA_seq = collapse_IgAseq(IgA_seq)
# DT:: datatable(IgA_seq@ig_seq_all, rownames = F) View(IgA_seq)


IgA_all = IgA_seq@ig_seq_all %>%
    separate(taxonomy, into = c("Reign", "Phylum", "Order", "Class", "Family", "Genus",
        "Species", "ASV", "rest"), sep = "#") %>%
    separate(col = new, into = colnames(sample_data(igaseq)), sep = "#")

IgA_all$alpha = ifelse(IgA_all$score > 1.96, "positively significative", ifelse(IgA_all$score <
    -1.96, "negatively significative", "not significative"))

dim(IgA_all)
#> [1] 4222   29
```

</div>

## Analyse and plot IgASeq

### Plot like Gordon IgASeq

<div style="text-align: justify">

For the example we will plot the IgAseq data like Gordon’s paper. We
first make a wilcoxon test to decipher which genera are significantly
different from zero. Then we plot as balloon plot using ggplot.

</div>

``` r
library(rstatix)
tmp = IgA_all %>%
    group_by(Genus, donor) %>%
    mutate(n = n()) %>%
    filter(n > 5, Genus != "NA") %>%
    wilcox_test(score ~ 1, mu = 0, alternative = "two.sided", detailed = T)

tmp2 = IgA_all %>%
    group_by(Genus, donor) %>%
    mutate(n = n()) %>%
    filter(n > 5, Genus != "NA") %>%
    dplyr::select(Phylum, Genus, score, sample_origin, donor) %>%
    mutate(mean_score = median(score)) %>%
    left_join(tmp)

tmp2$mean_score[tmp2$mean_score <= (-5)] = (-5)
alpha = ifelse(tmp2$p < 0.05, -log10(tmp2$p), 0.5)

tmp2 = tmp2 %>%
    ungroup() %>%
    select(p, score, mean_score, Phylum, Genus, donor) %>%
    mutate(alpha = ifelse(tmp2$p < 0.05, -log10(tmp2$p), 0), alpha2 = ifelse(tmp2$mean_score <
        0, -alpha, alpha), mean_score2 = scales::rescale(c(abs(mean_score)), to = c(0,
        5)))

p1 = tmp2 %>%
    ggplot(aes(0, Genus, fill = alpha2, size = mean_score2)) + geom_point(shape = 21) +
    scale_fill_gradient2(low = "olivedrab4", mid = "white", high = "brown", midpoint = 0,
        guide = F) + scale_size_continuous(breaks = c(0, 3, 3, 4, 4, 5, 5), labels = c(-5,
    -4, -3, 0, 3, 4, 5), range = c(0, 5)) + guides(size = guide_legend(override.aes = list(fill = colorRampPalette(c("darkgreen",
    "white", "brown"))(7), size = c(5, 4, 3, 1, 3, 4, 5)), nrow = 1, direction = "horizontal",
    title.position = "top", label.position = "bottom", label.hjust = 0.5, label.vjust = 1)) +
    labs(fill = "", x = "Age", size = "IgAseq median score") + facet_grid(Phylum ~
    donor, scales = "free", space = "free") + theme(axis.text.y = element_text(size = 8,
    face = "bold"), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    strip.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
    legend.position = "bottom", legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"), legend.box = "vertical",
    legend.box.background = element_rect(colour = "black"), strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"), panel.grid.major.y = element_line(linetype = 2,
        size = 0.05))
#> Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
#> ℹ Please use the `linewidth` argument instead.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.




p2 = tmp2 %>%
    ggplot(aes(x = score, Genus, fill = alpha2)) + geom_boxplot(outlier.size = 0) +
    scale_fill_gradient2(low = "olivedrab4", mid = "white", high = "brown", midpoint = 0,
        guide = F) + facet_grid(Phylum ~ donor, scales = "free_y", space = "free") +
    geom_vline(xintercept = 0, linetype = 2) + theme(axis.text.y = element_blank(),
    axis.text.x = element_blank(), strip.text.y = element_text(angle = 0, size = 10,
        face = "bold", hjust = 0), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(), legend.position = "bottom", legend.title = element_text(size = 10,
        face = "bold"), legend.text = element_text(size = 10, face = "bold"), legend.box = "vertical",
    legend.box.background = element_rect(colour = "black"), strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"), panel.grid.major.y = element_line(linetype = 2,
        size = 0.05))


ggpubr::ggarrange(p1, p2, nrow = 1, widths = c(1, 1), common.legend = T, legend = "bottom")
#> Warning: The `guide` argument in `scale_*()` cannot be `FALSE`. This was deprecated in
#> ggplot2 3.3.4.
#> ℹ Please use "none" instead.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
```

<div class="figure" style="text-align: center">

<img src="ImmuMicrobiome_book_files/figure-gfm/unnamed-chunk-22-1.png" alt="Bubble plot of IgASeq score"  />
<p class="caption">
Bubble plot of IgASeq score
</p>

</div>

### Reduction of dimensions on IgASeq

IgASeq is a particular type of data. Indeed, it is not compositional
anymore and only values outside \[-1.96 ; +1.96\] are considered as
statistically significant. What should be done for the taxa that are
comprised in this interval ? Good question.

No matter what, you cannot use classical distance algorithm such as
Bray-Curtis or Jaccard : those algorithm doesn’t handle negatie values.
We only have Euclidean and Canberra algorithm that can handle our new
type of data. It seems that Canberra handles well 0 centered values. So
we encourage the use of Canberra distance. The reduction approach is up
to you.

``` r
tmp = IgA_all %>%
    group_by(sample_origin, ASV, score, donor) %>%
    summarise(score = mean(score)) %>%
    pivot_wider(values_from = score, names_from = ASV, values_fill = 0) %>%
    as.data.frame

tmp %>%
    plot_reduction(method = "PCA", clinical_data = 1:2, dist = "canberra", group = "donor",
        type = "pure", draw = "polygon")

tmp %>%
    plot_reduction(method = "PCoA", clinical_data = 1:2, dist = "canberra", group = "donor",
        type = "pure", draw = "polygon")
```

<div class="figure" style="text-align: center">

<img src="ImmuMicrobiome_book_files/figure-gfm/unnamed-chunk-23-1.png" alt="Reduction of IgASeq score"  /><img src="ImmuMicrobiome_book_files/figure-gfm/unnamed-chunk-23-2.png" alt="Reduction of IgASeq score"  />
<p class="caption">
Reduction of IgASeq score
</p>

</div>

As you can see PCA can be strongly impacted by extreme values.
