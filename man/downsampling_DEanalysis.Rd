% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downsampling_DEanalysis.r
\name{downsampling_DEanalysis}
\alias{downsampling_DEanalysis}
\title{Downsample the dataset, based either on the individuals or cells, and run DE analysis on each downsampled output. Save results in a dataframe}
\usage{
downsampling_DEanalysis(
  SCE,
  range_downsampled = "placeholder",
  output_path = getwd(),
  sampled = "individuals",
  sampleID = "donor_id",
  design = "placeholder",
  sexID = "sex",
  celltypeID = "cell_type",
  assay_name = "counts",
  coef = "male",
  fdr = 0.05,
  nom_pval = 0.05,
  Nperms = 20,
  y = NULL,
  region = "single_region",
  control = NULL,
  pval_adjust_method = "BH",
  rmv_zero_count_genes = TRUE
)
}
\arguments{
\item{SCE}{the input data (should be an SCE object)}

\item{range_downsampled}{vector or list containing values which the SCE will be downsampled at, in ascending order}

\item{output_path}{base path in which outputs will be stored}

\item{sampled}{downsampling carried out based on what (either "individuals" or "cells")}

\item{sampleID}{sample ID}

\item{design}{the design formula of class type \code{formula}. Equation used to fit the model- data for the generalised linear model e.g. expression ~ sex + pmi + disease}

\item{sexID}{sex ID}

\item{celltypeID}{cell type ID}

\item{assay_name}{the assay name in the SCE object to use for the analysis. Default is "counts" which uses the counts assay in each SCE}

\item{coef}{which coefficient to carry out DE analysis with respect to}

\item{fdr}{the cut-off False Discovery Rate below which to select DEGs}

\item{nom_pval}{the cut-off nominal P-value below which to select DEGs (as an alternative to FDR)}

\item{Nperms}{number of subsets created when downsampling at each level}

\item{y}{the column name in the SCE object for the return variable e.g. "diagnosis" - Case or disease. Default is the last variable in the design formula. y can be discrete (logistic regression) or continuous (linear regression)}

\item{region}{the column name in the SCE object for the study region. If there are multiple regions in the study (for example two brain regions). Pseudobulk values can be derived separately. Default is "single_region" which will not split by region.}

\item{control}{character specifying which control level for the differential expression analysis e.g. in a case/control/other study use "control" in the y column to compare against. NOTE only need to specify if more than two groups in y, leave as default value for two groups or continuous y. Default is NULL.}

\item{pval_adjust_method}{the adjustment method for the p-value in the differential expression analysis. Default is benjamini hochberg "BH". See  stats::p.adjust for available options}

\item{rmv_zero_count_genes}{whether genes with no count values in any cell should be removed. Default is TRUE
Saves all DGE analysis outputs for downsampled files as well as a summary table of results showing number of true DEGs detected at each number of samples/cells}
}
\description{
Downsample the dataset, based either on the individuals or cells, and run DE analysis on each downsampled output. Save results in a dataframe
}
