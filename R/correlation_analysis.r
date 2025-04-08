#' Runs correlation analysis pipeline

#' @param dataset_name name of the dataset used to select significant DEGs from (specified as a string, name as in DEouts)
#' @param DEouts a list containing outputs of DGE analysis (as returned/optionally saved by DGE_analysis) for datasets to be used in the correlation analysis
#' @param celltype_correspondence list of different names specifying each cell type
#' @param pvals list of p-value cut-offs which will be used to select DEGs
#' @param data_names names of the datasets as they appear in the correlation plot
#' @param alphaval (alpha) transparency of the non-mean boxplots
#' @param N_randperms number of random permutations of the dataset used to select significant DEGs from
#' @param N_subsets number of pairs of random subsets of the dataset used to select significant DEGs from
#' @param sex_DEGs If TRUE, only keep genes present on sex chromosmomes. Queries hspanies gene Ensembl dataset.
#' @param fontsize_yaxislabels font size for axis labels in plot
#' @param fontsize_yaxisticks font size for axis tick labels in plot
#' @param fontsize_title font size for plot title
#' @param fontsize_legendlabels font size for legend labels in plot
#' @param fontsize_legendtitle font size for legend title in plot
#' @param fontsize_facet_labels font size for facet labels
#' @param output_path base path in which outputs will be stored

#' Saves all plots and DGE analysis outputs in the appropriate directories
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
#'
#' # Runnable example with bundled data
#' # 1. Prepare DGE
#' micro_tsai <- system.file("extdata", "Tsai_Micro.qs", package="poweranalysis")
#' SCE_tsai <- qs::qread(micro_tsai)
#' DGE_tsai <- poweranalysis::DGE_analysis(
#'     SCE_tsai,
#'     design = ~ sex,
#'     celltypeID="cluster_celltype",
#'     sampleID = "sample_id",
#'     coef = "M"
#' )
#' celltype_names <- list("Microglia" = list("Micro"))
#'
#' # 2. Run correlation analysis
#' correlation_analysis(
#'     dataset_name = "tsai",
#'     DEouts = list(DGE_tsai),
#'     celltype_correspondence = celltype_names,
#'     data_names = list("tsai"),
#'     sex_DEGs = TRUE,    # Keep only genes on sex chromosomes
#'     output_path = output_path
#' )
#'}

correlation_analysis <- function(dataset_name,
                                 DEouts,
                                 celltype_correspondence,
                                 pvals=c(0.05,0.025,0.01,0.001,0.0001),
                                 data_names="placeholder",
                                 alphaval=0.25,
                                 N_randperms=5,
                                 N_subsets=5,
                                 sex_DEGs=FALSE,
                                 fontsize_yaxislabels=12,
                                 fontsize_yaxisticks=9,
                                 fontsize_title=14,
                                 fontsize_legendlabels=9,
                                 fontsize_legendtitle=9,
                                 fontsize_facet_labels=9,
                                 output_path=getwd()){

    # Filter sex_DEGs
    if (sex_DEGs) {
        DEouts <- lapply(DEouts, function(x) {
            x$celltype_all_genes <- sex_chromosome_DEGs(x$celltype_all_genes)

            # Only keep DEGs for genes present in filtered all_genes
            celltype_DEGs <- lapply(names(x$celltype_DEGs), function(y) {
                celltype_DEG <- x$celltype_DEGs[[y]][x$celltype_DEGs[[y]]$name %in% x$celltype_all_genes[[y]]$name,]
                return(celltype_DEG)
            })
            names(celltype_DEGs) <- names(x$celltype_DEGs)
            x$celltype_DEGs <- celltype_DEGs
            return(x)
        })
    }


    # run plot_mean_correlation for each p-value (saving outputs)
    mean_correlation_results <- plot_mean_correlation(dataset_name=dataset_name,
                                                      DEouts=DEouts,
                                                      celltype_correspondence=celltype_correspondence,
                                                      pvals=pvals,
                                                      data_names=data_names,
                                                      output_path=output_path)

    # run correlation_boxplots for each p-value (saving outputs)
    correlation_boxplots(mean_correlation_results,
                         num_real_datasets=length(DEouts),
                         pvals=pvals,
                         N_randperms=N_randperms,
                         N_subsets=N_subsets,
                         sex_DEGs=sex_DEGs,
                         fontsize_yaxislabels=fontsize_yaxislabels,
                         fontsize_yaxisticks=fontsize_yaxisticks,
                         fontsize_title=fontsize_title,
                         fontsize_legendlabels=fontsize_legendlabels,
                         fontsize_legendtitle=fontsize_legendtitle,
                         fontsize_facet_labels=fontsize_facet_labels,
                         output_path=output_path)

}
