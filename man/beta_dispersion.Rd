% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/beta_dispersion.R
\name{beta_dispersion}
\alias{beta_dispersion}
\title{Make an unconstrained analysis of beta diversity.}
\usage{
beta_dispersion(
  physeq,
  axis_x = 1,
  axis_y = 2,
  dist = "bray",
  nf = 5,
  method = "PCoA",
  group = NULL,
  color_vector = c("cyan4", "brown", "deepskyblue", "black", "red"),
  legend_title = NULL,
  title = NULL,
  lwd = 1,
  conf = 0.9,
  cex = 2,
  font = 2,
  pch = 20,
  draw = "lines",
  ylimits = "auto",
  xlimits = "auto",
  text = F,
  ncol = 1,
  species = F,
  x.intersp = 1,
  y.intersp = 0.5,
  permanova = F,
  where = "topleft",
  inset = 0.2,
  boxplot = T,
  pca = T,
  ...
)
}
\arguments{
\item{physeq}{A phyloseq object.}

\item{axis_x}{Which component for axis x. Default = 1}

\item{axis_y}{Which component for axis y. Default = 2}

\item{dist}{Either a distance mecthod given as string or a distance object.
Default= "bray".}

\item{nf}{Number of component to keep. Default= 5.}

\item{method}{Method use to plot the beta diversity. Currently this function
supports only "PCoA".}

\item{group}{String defining the groups you want to plot on the beta diversity.}

\item{color_vector}{Optionnal. A color vector you wish to use for plotting}

\item{legend_title}{Optionnal. A legend title.}

\item{lwd}{Line width}

\item{conf}{Interval of confidence for the ellipses. Default= 0.9.}

\item{cex}{Text size}

\item{font}{Font type}

\item{pch}{Shape type}

\item{draw}{Either draw lines or polygon. Default = "lines".}

\item{ylimits}{Limits for Y axis, if not specified they will be automatically set}

\item{xlimits}{Limits for X axis, if not specified they will be automatically set}

\item{text}{Should the labels be printed}

\item{ncol}{Number of columns for the legend. Default = 1.}

\item{species}{Not implemented yet}

\item{x.intersp}{character interspacing factor for horizontal (x) spacing between symbol and legend text.}

\item{y.intersp}{vertical (y) distances (in lines of text shared above/below each legend entry). A vector with one element for each row of the legend can be used.}

\item{permanova}{Make the permanova analysis. Default = FALSE.}

\item{where}{Where you want the permanova result.}

\item{inset}{Permanova legend adjustment see \link{legend}}

\item{boxplot}{Draw boxplot for each axis. Default = T}

\item{pca}{Argument for t-SNE, does the t-SNE needs to run on a PCA first or not ? Default=FALSE}
}
\description{
This a wrapper function to plot beta diversity directly from a phyloseq
object. This is based on the base package and mostly \link{vegan}.
}
