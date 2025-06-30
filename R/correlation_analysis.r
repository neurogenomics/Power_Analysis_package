#' Perform correlation analysis of DEG effect sizes between single-cell datasets
#'
#' Runs the correlation analysis pipeline by computing Spearman’s rank correlations of log₂ fold-changes for differentially expressed genes (DEGs) across and within multiple scRNA-seq datasets. Uses a user-specified reference dataset to define DEGs and compares effect sizes across studies, independently sampled subsets, and permuted controls.

#' @param main_dataset Name of the dataset used to select significant DEGs from (specified as a string, use the dataset name as in dataset_names)
#' @param SCEs A list of SingleCellExperiment (SCE) objects, each representing a scRNA-seq dataset.
#' @param sampleIDs A character vector specifying the column name in each SCE that represents sample or donor IDs (in order of SCEs).
#' @param celltypeIDs A character vector specifying the column name in each SCE that denotes cell type identity (in order of SCEs).
#' @param celltype_correspondence A named vector that maps a standard cell type label (e.g., list(Micro=c("Micro",NA), Astro=c(NA,"Astro")) to how that cell type appears in each dataset. Use `NA` if the cell type is not present in a given dataset.
#' @param dataset_names A vector of names corresponding to each dataset (as you would like them to appear in output plots).
#' @param assay_names A character vector specifying the assay names in each SCE that will be used for the analysis (in order of SCEs). Default is a vector with all entries `"counts"`, which uses the count assay in each SCE.
#' @param pvals list of P-value thresholds for selecting DEGs in each individual dataset. Default is c(0.05,0.025,0.01,0.001,0.0001).
#' @param alphaval Transparency of the non-mean boxplots. The value of alpha ranges between 0 (completely transparent) and 1 (completely opaque).
#' @param N_randperms Number of random permutations of the dataset used to select significant DEGs from. Default is 5.
#' @param N_subsets Number of pairs of random subsets of the dataset used to select significant DEGs from. Default is 5.
#' @param output_path A directory path where outputs will be saved.
#' @param sex_DEGs If TRUE, only keep genes present on sex chromosmomes. Queries hspanies gene Ensembl dataset.
#' @param fontsize_yaxislabels font size for axis labels in plot
#' @param fontsize_yaxisticks font size for axis tick labels in plot
#' @param fontsize_title font size for plot title
#' @param fontsize_legendlabels font size for legend labels in plot
#' @param fontsize_legendtitle font size for legend title in plot
#' @param fontsize_facet_labels font size for facet labels

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
#'                                "Astrocytes" = c("Astrocytes", "astro", NA),
#'                                "Oligodendrocytes" = c("oligo", "Oligodendrocytes", NA))
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
#'     main_dataset = "tsai",
#'     SCEs = list(tsai),
#'     celltype_correspondence = celltype_names,
#'     dataset_names = list("tsai"),
#'     sex_DEGs = TRUE,    # Keep only genes on sex chromosomes
#'     output_path = output_path
#' )
#'}

correlation_analysis <- function(main_dataset,
                                 SCEs,
                                 sampleIDs,
                                 celltypeIDs,
                                 celltype_correspondence,
                                 dataset_names,
                                 assay_names="counts",
                                 pvals=c(0.05,0.025,0.01,0.001,0.0001),
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

    # Comprehensive validation for all parameters used in the pipeline
    validate_input_parameters_correlation(main_dataset=main_dataset,
                                          SCEs=SCEs,
                                          sampleIDs=sampleIDs,
                                          celltypeIDs=celltypeIDs,
                                          celltype_correspondence=celltype_correspondence,
                                          dataset_names=dataset_names,
                                          assay_names=assay_names,
                                          pvalues=pvals,
                                          alphaval=alphaval,
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

    # run plot_mean_correlation for each p-value (saving outputs)
    mean_correlation_results <- plot_mean_correlation(main_dataset=main_dataset,
                                                      SCEs=SCEs,
                                                      sampleIDs=sampleIDs,
                                                      celltypeIDs=celltypeIDs,
                                                      celltype_correspondence=celltype_correspondence,
                                                      dataset_names=dataset_names,
                                                      assay_names=assay_names,
                                                      pvals=pvals,
                                                      N_randperms=N_randperms,
                                                      N_subsets=N_subsets,
                                                      sex_DEGs=sex_DEGs,
                                                      output_path=output_path)

    # run correlation_boxplots for each p-value (saving outputs), unless only one SCE
    if(length(SCEs) == 1){
        print("Only one SCE provided, skipping correlation boxplots")
    }else{
        correlation_boxplots(mean_correlation_results,
                             num_real_datasets=length(SCEs),
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

}
