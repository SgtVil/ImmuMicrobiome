% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_all_variance.R
\name{plot_all_variance}
\alias{plot_all_variance}
\title{Plot the overall variance of your dataset}
\usage{
plot_all_variance(
  variance,
  plot_type = "boxplot",
  top = 30,
  col = c("brown", "orange", "grey")
)
}
\arguments{
\item{variance}{An object returned by the \link{calculate_variance}}

\item{plot_type}{Categorical. Plot as boxplots or heatmap. Default= "boxplot".}

\item{top}{For heatmap plotting only. The top N features to plot. Default = 30.}
}
\value{
A ggplot.
}
\description{
This function will plot the overall variance explained by the selected factors.
You can choose to plot it as boxplots or heatmap.
\itemize{
\item Boxplots : points will represent the mean abundance of the features. Features with p.adj < 0.05 and p < 0.05 will be plotted.
\item Heatmap : all factors will be plotted as rows. Currently doesn't shows which features ar significant.
}
}
\examples{
var = metabolomic \%>\%
dplyr::select(!child_id) \%>\%
calculate_variance(clinical_data = 1:3, cores = 1)

# A boxplot with jitter.
plot_all_variance(var, col = c("brown", "darkgreen", "grey"))

# alternatively you can plot the results as a Heatmap
plot_all_variance(var, plot_type = "heatmap")
}
