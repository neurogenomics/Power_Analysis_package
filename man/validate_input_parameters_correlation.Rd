% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validate_input_parameters_correlation.r
\name{validate_input_parameters_correlation}
\alias{validate_input_parameters_correlation}
\title{Tests input parameters for functions}
\usage{
validate_input_parameters_correlation(
  dataset_name = "placeholder",
  allstudies = "placeholder",
  celltypes = "placeholder",
  pvalue = "placeholder",
  data_names = "placeholder",
  corrMats = "placeholder",
  numRealDatasets = "placeholder",
  alphaval = "placeholder",
  numPerms = "placeholder",
  numSubsets = "placeholder",
  sexDEGs = "placeholder",
  fontsize_yaxislabels = "placeholder",
  fontsize_yaxisticks = "placeholder",
  fontsize_title = "placeholder",
  fontsize_legendlabels = "placeholder",
  fontsize_legendtitle = "placeholder",
  fontsize_facet_labels = "placeholder",
  output_path = "placeholder"
)
}
\arguments{
\item{dataset_name}{name of the dataset used to select significant DEGs from (specified as a string, name as in allStudies)}

\item{allstudies}{a list containing all the datasets (most likely as SCE objects)}

\item{celltypes}{a list containing the celltypes to compute mean correlation across}

\item{pvalue}{the cut-off p-value which will be used to select DEGs}

\item{data_names}{names of the datasets as they appear in the correlation plot}

\item{corrMats}{(named) list of correlation matrices for each celltype with the final element being the mean correlation matrix, all at specified p-value}

\item{numRealDatasets}{total number of \emph{real} datasets (most likely the number of studies, but sometimes a study may be split e.g. into 2 brain regions, so in this case it would be the number of studies plus 1)}

\item{alphaval}{(alpha) transparency of the non-mean boxplots}

\item{numPerms}{number of random permutations of the dataset used to select significant DEGs from}

\item{numSubsets}{number of pairs of random subsets of the dataset used to select significant DEGs from}

\item{sexDEGs}{true if DEGs come from sex chromosomes, else false}

\item{fontsize_yaxislabels}{font size for axis labels in plot}

\item{fontsize_yaxisticks}{font size for axis tick labels in plot}

\item{fontsize_title}{font size for plot title}

\item{fontsize_legendlabels}{font size for legend labels in plot}

\item{fontsize_legendtitle}{font size for legend title in plot}

\item{fontsize_facet_labels}{font size for facet labels}

\item{output_path}{base path in which outputs will be stored
Checks all correlation analysis parameters are specified correctly}
}
\description{
Tests input parameters for functions
}
