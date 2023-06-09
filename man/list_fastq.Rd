% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/list_fastq.R
\name{list_fastq}
\alias{list_fastq}
\title{Make the list of fastq files for the rest of the proceedings.}
\usage{
list_fastq(path, pattern, separator = "_", level)
}
\arguments{
\item{path}{Path of selected samples}

\item{pattern}{Pattern for forward and reverse fastq if paired sample analysis, single pattern for single end.}

\item{separator}{Separator by which you want to separate the filenames to obtain sample names. Default = "_".}

\item{level}{Which level of the separated filenames you want to keep as sample names.}
}
\value{
A list of fwd, rv and names.
}
\description{
This will enable other wrapper functions for dada2 pipeline.
}
\examples{
No example
}
