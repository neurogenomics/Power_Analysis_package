#' Perform cross-dataset power analysis using scRNA-seq and bulk RNA-seq DEG overlap
#'
#' Runs the complete bulk RNA-seq-informed power analysis pipeline by performing downsampling-based DEG detection across multiple scRNA-seq datasets, comparing overlaps with bulk RNA-seq DEGs, and generating summary plots for evaluation.

#' @param SCEs A list of SingleCellExperiment (SCE) objects, each representing a scRNA-seq dataset.
#' @param dataset_names A vector of names corresponding to each dataset (as you would like them to appear in output plots).
#' @param celltype_correspondence A named vector that maps a standard cell type label (e.g., list(Micro=c("Micro",NA), Astro=c(NA,"Astro")) to how that cell type appears in each dataset. Use `NA` if the cell type is not present in a given dataset.
#' @param output_path A clean directory path where DGE analysis outputs of down-sampled datasets and summary plots will be saved (should contain no subdirectories).
#' @param celltypeIDs A character vector specifying the column name in each SCE that denotes cell type identity (in order of SCEs).
#' @param sampleIDs  A character vector specifying the column name in each SCE that represents sample or donor IDs (in order of SCEs).
#' @param sampled Specifies the unit of down-sampling. Can be either `"individuals"` or `"cells"`, depending on whether the analysis downsamples across samples or cells.
#' @param bulkDE DGE analysis output for a bulk RNA-seq dataset (e.g., `LFSR.tsv`): rows (rownames) should be the genes, columns should be tissues, and entries should be significance levels
#' @param bulk_cutoff Proportion (0–1) of bulk tissues in which a gene must be differentially expressed to be considered (e.g., 0.9 selects DEGs found in ≥90% of tissues). Default is 0.9.
#' @param pvalue P-value threshold for selecting DEGs in each individual dataset. Default is 0.05.
#' @param Nperms Number of subsets (permutations) to generate at each downsampling level during power analysis. Each subset is analyzed independently to estimate variability. Default is 20.
#' @param fontsize_axislabels Font size for axis labels in plot
#' @param fontsize_axisticks Font size for axis tick labels in plot
#' @param fontsize_title Font size for plot title
#' @param fontsize_legendlabels Font size for legend labels in plot
#' @param fontsize_legendtitle Font size for legend title in plot
#' @param plot_title Plot title

#' Saves all plots in the appropriate directories
#' @export

#' @examples
#'\dontrun{
#'# Example of celltype_correspondence:
#'celltype_correspondence <- list([overall_celltype_name1] = c([celltype1.1], [celltype1.2], ...),
#'                                 [overall_celltype_name2] = c([celltype2.1], [celltype2.2], ...),
#'                                 ...)
#'# E.g.
#'celltype_correspondence <- list("Microglia" = c("micro", "microglial cells", "Micro"),
#'                                "Astrocytes" = c("Astrocytes", "astro"),
#'                                "Oligodendrocytes" = c("oligo", "Oligodendrocytes"))
#'# If a celltype is not present in the dataset, the corresponding value should be set to NA. E.g.:
#'celltype_correspondence <- list("Microglia" = c("micro",NA), "Endothelial" = c(NA,"endo"))
#'}

bulk_power_analysis <- function(SCEs,
                                dataset_names,
                                celltype_correspondence,
                                output_path=getwd(),
                                celltypeIDs="cell_type",
                                sampled="individuals",
                                sampleIDs="donor_id",
                                bulkDE="placeholder",
                                bulk_cutoff=0.9,
                                pvalue=0.05,
                                Nperms=20,
                                fontsize_axislabels=12,
                                fontsize_axisticks=9,
                                fontsize_title=14,
                                fontsize_legendlabels=9,
                                fontsize_legendtitle=9,
                                plot_title="placeholder"){

    # Run bulk_downsampling_DGEanalysis for all cell types
    bulk_downsampling_DGEanalysis(SCEs = SCEs,
                                  dataset_names = dataset_names,
                                  celltype_correspondence = celltype_correspondence,
                                  sampled = sampled,
                                  sampleIDs = sampleIDs,
                                  celltypeIDs = celltypeIDs,
                                  output_path = output_path,
                                  pvalue = pvalue,
                                  Nperms = Nperms)

    # Run gather_celltype_DEGs
    gather_celltype_DEGs(celltype_correspondence = celltype_correspondence,
                         pvalue = pvalue,
                         Nperms = Nperms,
                         output_path = output_path,
                         sampled = sampled)


    # Run prop_bulk_DEGs_sc
    prop_bulk_DEGs_sc(bulkDE = bulkDE,
                      output_path = output_path,
                      bulk_cutoff = bulk_cutoff,
                      pvalue = pvalue,
                      sampled = sampled)

}
