% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differential_abundance.R
\name{differential_abundance}
\alias{differential_abundance}
\title{Title
Find differential abundance within a dataset.}
\usage{
differential_abundance(
  physeq,
  group,
  paired = F,
  col1 = "Genus",
  col2 = NULL,
  cores = F,
  plot = F
)
}
\arguments{
\item{physeq}{Phyloseq object to compute analysis.}

\item{group}{The group to perform analysis. Make sure the vector is of length 2. If not please subset your data.}

\item{paired}{Are the data paired.}

\item{col1}{First column used to name the taxa in the plot. Default is Genus level.}

\item{col2}{Second column used to name the taxa in the plot. Default is the last column of your
tax_table. It is recommended to use an ASV or OTU number for each taxa.}
}
\value{
A list of three elements.
\item{all_features}{A dataframe containing all taxa tested with their respective values.}
\item{signif_features}{A dataframe containing only taxa significant for eBH test.}
\item{barplot}{Differential abundance plotted witha barplot}
\item{volcano}{Differential abundance plotted witha volcano plot}
}
\description{
This is a wrapper function that use \link{ALDEx2} package to find the taxa that are differentially abundant
between two groups. You have to provide a phyloseq object with only two groups in order to use this function.
}
