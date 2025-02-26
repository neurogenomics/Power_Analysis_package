% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validate_input_parameters_power.r
\name{validate_input_parameters_power}
\alias{validate_input_parameters_power}
\title{Tests input parameters for functions}
\usage{
validate_input_parameters_power(
  data = "placeholder",
  range_downsampled = "placeholder",
  output_path = "placeholder",
  inpath = "placeholder",
  sampled = "placeholder",
  sampleID = "placeholder",
  design = "placeholder",
  sexID = "placeholder",
  celltypeID = "placeholder",
  coeff = "placeholder",
  fdr = "placeholder",
  nom_pval = "placeholder",
  Nperms = "placeholder",
  N_randperms = "placeholder",
  N_subsetpairs = "placeholder",
  y = "placeholder",
  region = "placeholder",
  control = "placeholder",
  pval_adjust_method = "placeholder",
  rmv_zero_count_genes = "placeholder",
  base_outpath = "placeholder",
  datasets_path = "placeholder",
  sampleIDs = "placeholder"
)
}
\arguments{
\item{data}{the input data (should be an SCE object)}

\item{range_downsampled}{vector or list containing values which the data will be downsampled at, in ascending order}

\item{output_path}{base path in which outputs will be stored}

\item{inpath}{base path where downsampled DGE analysis output is stored (taken to be output_path if not provided)}

\item{sampled}{downsampling carried out based on what (either "individuals" or "cells")}

\item{sampleID}{sample ID}

\item{design}{the design formula of class type \code{formula}. Equation used to fit the model- data for the generalised linear model e.g. expression ~ sex + pmi + disease}

\item{sexID}{sex ID}

\item{celltypeID}{cell type ID}

\item{coeff}{which coefficient to carry out DE analysis with respect to}

\item{fdr}{the cut-off False Discovery Rate below which to select DEGs}

\item{nom_pval}{the cut-off nominal P-value below which to select DEGs (as an alternative to FDR)}

\item{Nperms}{number of subsets created when downsampling at each level}

\item{N_randperms}{number of randomised permutations of the dataset (based on sex) to be correlated}

\item{N_subsetpairs}{number of pairs of subsets of the dataset to be correlated}

\item{y}{the column name in the SCE object for the return variable e.g. "diagnosis" - Case or disease. Default is the last variable in the design formula. y can be discrete (logistic regression) or continuous (linear regression)}

\item{region}{the column name in the SCE object for the study region. If there are multiple regions in the study (for example two brain regions). Pseudobulk values can be derived separately. Default is "single_region" which will not split by region.}

\item{control}{character specifying which control level for the differential expression analysis e.g. in a case/control/other study use "control" in the y column to compare against. NOTE only need to specify if more than two groups in y, leave as default value for two groups or continuous y. Default is NULL.}

\item{pval_adjust_method}{the adjustment method for the p-value in the differential expression analysis. Default is benjamini hochberg "BH". See  stats::p.adjust for available options}

\item{rmv_zero_count_genes}{whether genes with no count values in any cell should be removed. Default is TRUE}

\item{base_outpath}{the base output path which contains power analysis output folders for \emph{all} datasets (each one as created by "downsampling_DEanalysis.r")}

\item{datasets_path}{path where all dataset SCEs (in .qs format) are saved}

\item{sampleIDs}{list of sample ID names for datasets in "datasets_path" (may be in any order)
Checks all power analysis parameters are specified correctly}
}
\description{
Tests input parameters for functions
}
