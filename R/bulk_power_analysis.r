#' Runs entire bulk RNA-seq power analysis pipeline

#' @param SCEs list of the input data (elements should be SCE objects)
#' @param celltype the cell type we are focusing on (name as it appears in cell type sub-directory name)
#' @param celltype_correspondence list of different names specifying each cell type
#' @param output_path path storing the down-sampled DGE analysis outputs for each dataset (and saves outputs)
#' @param range_downsampled vector or list containing values which the data will be downsampled at, in ascending order
#' @param celltypeID cell type ID
#' @param sampled downsampling carried out based on what (either "individuals" or "cells")
#' @param sampleID sample ID
#' @param bulkDE DGE analysis output for a bulk RNA-seq dataset: rows (rownames) should be the genes, columns should be tissues, and entries should be significance levels
#' @param bulk_cutoff percentage (proportion), specified so that we select DEGs common across >=bulk_cutoff of the tissues in the Bulk dataset
#' @param pvalue the cut-off p-value used to select DEGs (for bulk data)
#' @param Nperms number of permutations of DGE analysis outputs for each sample
#' @param fontsize_axislabels font size for axis labels in plot
#' @param fontsize_axisticks font size for axis tick labels in plot
#' @param fontsize_title font size for plot title
#' @param fontsize_legendlabels font size for legend labels in plot
#' @param fontsize_legendtitle font size for legend title in plot
#' @param plot_title plot title

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
#'}

bulk_power_analysis <- function(SCEs,
                                celltype,
                                celltype_correspondence,
                                output_path=getwd(),
                                range_downsampled="placeholder",
                                celltypeID="cell_type",
                                sampled="individuals",
                                sampleID="donor_id",
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
                                  sampled = sampled,
                                  sampleID = sampleID,
                                  celltypeID = celltypeID,
                                  celltype_correspondence = celltype_correspondence,
                                  output_path = output_path)

    # Run gather_celltype_DEGs
    gather_celltype_DEGs(output_path = output_path,
                         range_downsampled = range_downsampled,
                         celltype_correspondence = celltype_correspondence,
                         Nperms = Nperms,
                         pvalue = pvalue)

    # Run prop_bulk_DEGs_sc
    prop_bulk_DEGs_sc(bulkDE = bulkDE,
                      output_path = output_path,
                      range_downsampled = range_downsampled,
                      bulk_cutoff = bulk_cutoff,
                      pvalue = pvalue)

    # Run prop_bulk_DEGs_sc_celltype
    prop_bulk_DEGs_sc_celltype(bulkDE = bulkDE,
                               output_path = output_path,
                               range_downsampled = range_downsampled,
                               celltype = celltype,
                               celltype_correspondence = celltype_correspondence,
                               pvalue = pvalue)

}