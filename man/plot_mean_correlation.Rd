% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_mean_correlation.r
\name{plot_mean_correlation}
\alias{plot_mean_correlation}
\title{Obtain the average correlation (across celltypes) at a specified cutoff p-value}
\usage{
plot_mean_correlation(
  dataset_name,
  allstudies,
  celltypes,
  pvalue,
  data_names = "placeholder"
)
}
\arguments{
\item{dataset_name}{name of the dataset used to select significant DEGs from (specified as a string, name as in allStudies)}

\item{allstudies}{a list containing all the datasets (most likely as SCE objects)}

\item{celltypes}{a list containing the celltypes to compute mean correlation across}

\item{pvalue}{the cut-off p-value which will be used to select DEGs}

\item{data_names}{names of the datasets as they appear in the correlation plot}
}
\value{
mean correlation matrix
}
\description{
Obtain the average correlation (across celltypes) at a specified cutoff p-value
}
