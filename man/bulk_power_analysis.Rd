% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bulk_power_analysis.r
\name{bulk_power_analysis}
\alias{bulk_power_analysis}
\title{Perform cross-dataset power analysis using scRNA-seq and bulk RNA-seq DEG overlap}
\usage{
bulk_power_analysis(
  SCEs,
  dataset_names,
  celltype_correspondence,
  output_path = getwd(),
  celltypeIDs = "cell_type",
  sampled = "individuals",
  sampleIDs = "donor_id",
  assay_names = "counts",
  bulkDE = "placeholder",
  bulk_cutoff = 0.9,
  pvalue = 0.05,
  Nperms = 20,
  fontsize_axislabels = 12,
  fontsize_axisticks = 9,
  fontsize_title = 14,
  fontsize_legendlabels = 9,
  fontsize_legendtitle = 9,
  plot_title = "placeholder"
)
}
\arguments{
\item{SCEs}{A list of SingleCellExperiment (SCE) objects, each representing a scRNA-seq dataset.}

\item{dataset_names}{A vector of names corresponding to each dataset (as you would like them to appear in output plots).}

\item{celltype_correspondence}{A named vector that maps a standard cell type label (e.g., list(Micro=c("Micro",NA), Astro=c(NA,"Astro")) to how that cell type appears in each dataset. Use \code{NA} if the cell type is not present in a given dataset.}

\item{output_path}{A clean directory path where DGE analysis outputs of down-sampled datasets and summary plots will be saved (should contain no subdirectories).}

\item{celltypeIDs}{A character vector specifying the column name in each SCE that denotes cell type identity (in order of SCEs).}

\item{sampled}{Specifies the unit of down-sampling. Can be either \code{"individuals"} or \code{"cells"}, depending on whether the analysis downsamples across samples or cells.}

\item{sampleIDs}{A character vector specifying the column name in each SCE that represents sample or donor IDs (in order of SCEs).}

\item{assay_names}{A character vector specifying the assay names in each SCE that will be used for the analysis (in order of SCEs). Default is a vector with all entries \code{"counts"}, which uses the count assay in each SCE.}

\item{bulkDE}{DGE analysis output for a bulk RNA-seq dataset (e.g., \code{LFSR.tsv}): rows (rownames) should be the genes, columns should be tissues, and entries should be significance levels}

\item{bulk_cutoff}{Proportion (0–1) of bulk tissues in which a gene must be differentially expressed to be considered (e.g., 0.9 selects DEGs found in ≥90\% of tissues). Default is 0.9.}

\item{pvalue}{P-value threshold for selecting DEGs in each individual dataset. Default is 0.05.}

\item{Nperms}{Number of subsets (permutations) to generate at each downsampling level during power analysis. Each subset is analyzed independently to estimate variability. Default is 20.}

\item{fontsize_axislabels}{Font size for axis labels in plot}

\item{fontsize_axisticks}{Font size for axis tick labels in plot}

\item{fontsize_title}{Font size for plot title}

\item{fontsize_legendlabels}{Font size for legend labels in plot}

\item{fontsize_legendtitle}{Font size for legend title in plot}

\item{plot_title}{Plot title
Saves all plots in the appropriate directories}
}
\description{
Runs the complete bulk RNA-seq-informed power analysis pipeline by performing downsampling-based DEG detection across multiple scRNA-seq datasets, comparing overlaps with bulk RNA-seq DEGs, and generating summary plots for evaluation.
}
\examples{
\dontrun{
# Example of celltype_correspondence:
celltype_correspondence <- list([overall_celltype_name1] = c([celltype1.1], [celltype1.2], ...),
                                [overall_celltype_name2] = c([celltype2.1], [celltype2.2], ...),
                                ...)
# E.g.
celltype_correspondence <- list("Microglia" = c("micro", "microglial cells", "Micro"),
                               "Astrocytes" = c("Astrocytes", "astro"),
                               "Oligodendrocytes" = c("oligo", "Oligodendrocytes"))
# If a celltype is not present in the dataset, the corresponding value should be set to NA. E.g.:
celltype_correspondence <- list("Microglia" = c("micro",NA), "Endothelial" = c(NA,"endo"))
}
}
