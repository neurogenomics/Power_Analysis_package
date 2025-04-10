% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prop_bulk_DEGs_sc_celltype.r
\name{prop_bulk_DEGs_sc_celltype}
\alias{prop_bulk_DEGs_sc_celltype}
\title{Obtain percentage overlap between DEGs from bulk data (DE across all tissues) and various scRNA-seq datasets, for a specified cell type}
\usage{
prop_bulk_DEGs_sc_celltype(
  bulkDE,
  range_downsampled,
  celltype_correspondence,
  celltype,
  sampled = "individuals",
  bulk_cutoff = 0.9,
  pvalue = 0.05,
  fontsize_axislabels = 12,
  fontsize_axisticks = 9,
  fontsize_title = 14,
  fontsize_legendlabels = 9,
  fontsize_legendtitle = 9,
  plot_title = "placeholder",
  output_path = getwd()
)
}
\arguments{
\item{bulkDE}{DGE analysis output for a bulk RNA-seq dataset: rows (rownames) should be the genes, columns should be tissues, and entries should be significance levels}

\item{range_downsampled}{vector or list containing values which the data will be downsampled at, in ascending order}

\item{celltype_correspondence}{list of different names specifying each cell type}

\item{celltype}{the cell type we are focusing on (name as it appears in cell type sub-directory name)}

\item{sampled}{downsampling carried out based on what (either "individuals" or "cells")}

\item{bulk_cutoff}{percentage (proportion between 0 and 1), specified so that we select DEGs common across >= bulk_cutoff of the tissues in the Bulk dataset}

\item{pvalue}{the cut-off p-value used to select DEGs (for both, bulk and scRNA-seq datasets)}

\item{fontsize_axislabels}{font size for axis labels in plot}

\item{fontsize_axisticks}{font size for axis tick labels in plot}

\item{fontsize_title}{font size for plot title}

\item{fontsize_legendlabels}{font size for legend labels in plot}

\item{fontsize_legendtitle}{font size for legend title in plot}

\item{plot_title}{plot title}

\item{output_path}{path storing the down-sampled DGE analysis for each single-cell dataset
Saves plot showing percentage DEGs from bulk data found in each scRNA-seq dataset, for a specified cell type, in the appropriate directory}
}
\description{
Obtain percentage overlap between DEGs from bulk data (DE across all tissues) and various scRNA-seq datasets, for a specified cell type
}
