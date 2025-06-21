# Define global variables
utils::globalVariables(c("dataset"))


#' Perform robust power analysis for differential gene expression in scRNA-seq dataset
#'
#' Run the complete power analysis pipeline by downsampling individuals and cells, performing differential expression analysis, and generating power plots.

#' @importFrom stats as.formula

#' @param SCE A `SingleCellExperiment` object containing the input scRNA-seq data. You may also provide a path to an `.R`, `.rds`, or `.qs` file. If using a file, ensure the `SCE` object inside is named `SCE`.
#' @param range_downsampled_individuals A numeric vector specifying the number of individuals to include at each downsampling level, in ascending order (e.g., c(10, 20, 30)). By default, 12 evenly spaced values are generated from 0 to the total number of samples, each rounded up to the nearest multiple of 5.
#' @param range_downsampled_cells A numeric vector specifying the number of cells per individual to include at each downsampling level, in ascending order (e.g., c(20, 40, 60)). By default, 11 evenly spaced values are generated from 0 to the 90th percentile of per-individual cell counts, each rounded to the nearest multiple of 5.
#' @param output_path A directory path where DGE analysis outputs of down-sampled datasets and power plots will be saved.
#' @param sampleID Name of the column in the `SCE` metadata that identifies biological replicates (e.g., patient ID). This column is used for grouping in the pseudobulk approach.
#' @param celltypeID Name of the column in the `SCE` metadata indicating cell type labels. This is used to identify celltype specific DEGs.
#' @param sexID Name of the column in the `SCE` metadata that encodes the sex of individuals. Default is `"sex"`.
#' @param design  A model formula specifying covariates for differential expression analysis. It should be of class `formula` (e.g., `~ sex + pmi + disease`). This formula is used to fit a generalized linear model.
#' @param y Name of the column in the `SCE` metadata representing the response variable (e.g., "diagnosis" - case or disease). If not specified, defaults to the last variable in the `design` formula. Accepts both categorical (logistic regression) and continuous (linear regression) variables.
#' @param coeff Character string indicating the level of the response variable (`y`) to test for in differential expression. For case-control studies, this would typically be "case" (e.g. "AD"). Typically used in binary comparisons. Not required for continuous outcomes.
#' @param fdr Adjusted p-value (False Discovery Rate) threshold for selecting significantly differentially expressed genes (DEGs). Only genes with adjusted p-values below this value will be retained. Default is 0.05.
#' @param nom_pval Nominal (unadjusted) p-value threshold for selecting DEGs. Used as an alternative to FDR when preferred. Only genes with p-values below this cutoff will be retained. Default is 0.05.
#' @param Nperms Number of subsets (permutations) to generate at each downsampling level during power analysis. Each subset is analyzed independently to estimate variability. Default is 20.
#' @param region Optional column in `SCE` metadata indicating the tissue or brain region. If present, differential expression is performed within each region separately. Defaults to "single_region" (i.e., no regional split).
#' @param control  Optional. Character string specifying the control level in the response variable (`y`) to compare against. Only required if `y` contains more than two levels. Ignored for binary or continuous outcomes.
#' @param pval_adjust_method Method used to adjust p-values for multiple testing. Default is "BH" (Benjamini–Hochberg). See `stats::p.adjust` for available options.
#' @param rmv_zero_count_genes Logical. Whether to remove genes with zero counts across all cells. Default is `TRUE`.

#' Saves all plots and DGE analysis outputs in the appropriate directories
#' @export
#'
#' @examples
#'\dontrun{
#' # Too slow to run with check()
#' # 1. Prepare SCE
#' micro_tsai <- system.file("extdata", "Tsai_Micro.qs", package="poweranalysis")
#' SCE_tsai <- qs::qread(micro_tsai)
#'
#' # 2. Run Power Analysis
#' PA_tsai <- poweranalysis::power_analysis(
#'     SCE_tsai,
#'     sampleID = "sample_id",
#'     celltypeID = "cluster_celltype",
#'     design = ~ sex,
#'     coef = "M",
#'     output_path = tempdir()
#' )
#' PA_tsai
#'}
#'
#'

power_analysis <- function(SCE,
                           range_downsampled_individuals="placeholder",
                           range_downsampled_cells="placeholder",
                           output_path=getwd(),
                           sampleID="donor_id",
                           design="placeholder",
                           sexID="sex",
                           celltypeID="cell_type",
                           coeff="male",
                           fdr=0.05,
                           nom_pval=0.05,
                           Nperms=20,
                           y=NULL,
                           region="single_region",
                           control=NULL,
                           pval_adjust_method="BH",
                           rmv_zero_count_genes=TRUE){

    # Comprehensive validation for all parameters used in the pipeline
    validate_input_parameters_power(SCE=SCE,
                                    range_downsampled=range_downsampled_individuals,
                                    output_path=output_path,
                                    sampleID=sampleID,
                                    design=design,
                                    sexID=sexID,
                                    celltypeID=celltypeID,
                                    coeff=coeff,
                                    fdr=fdr,
                                    nom_pval=nom_pval,
                                    Nperms=Nperms,
                                    y=y,
                                    region=region,
                                    control=control,
                                    pval_adjust_method=pval_adjust_method,
                                    rmv_zero_count_genes=rmv_zero_count_genes)

    setwd(output_path)

    # alter range_downsampled_individuals
    if(identical(range_downsampled_individuals,"placeholder")){
        range_downsampled_individuals <- downsampling_range(SCE, "individuals", sampleID)
    }
    # alter range_downsampled_cells
    if(identical(range_downsampled_cells,"placeholder")){
        range_downsampled_cells <- downsampling_range(SCE, "cells", sampleID)
    }
    # alter design
    if(design=="placeholder"){
        design=as.formula(paste0("~",sexID))
    }

    # create preliminary plots
    preliminary_plots(SCE=SCE,
                      output_path=output_path,
                      sampleID=sampleID,
                      design=design,
                      sexID=sexID,
                      celltypeID=celltypeID,
                      coeff=coeff,
                      fdr=fdr)

    ## create power plots for down-sampling individuals, cells
    # down-sample individuals and run DE analysis
    downsampling_DEanalysis(SCE=SCE,
                            range_downsampled=range_downsampled_individuals,
                            output_path=output_path,
                            sampled="individuals",
                            sampleID=sampleID,
                            design=design,
                            sexID=sexID,
                            celltypeID=celltypeID,
                            y=y,
                            region=region,
                            control=control,
                            pval_adjust_method=pval_adjust_method,
                            rmv_zero_count_genes=rmv_zero_count_genes,
                            coeff=coeff,
                            fdr=fdr,
                            nom_pval=nom_pval,
                            Nperms=Nperms)
    # create power plots
    power_plots(SCE=SCE,
                range_downsampled=range_downsampled_individuals,
                output_path=output_path,
                sampled="individuals",
                sampleID=sampleID,
                celltypeID=celltypeID,
                fdr=fdr,
                nom_pval=nom_pval,
                Nperms=Nperms)
    # down-sample cells and run DE analysis
    downsampling_DEanalysis(SCE=SCE,
                            range_downsampled=range_downsampled_cells,
                            output_path=output_path,
                            sampled="cells",
                            sampleID=sampleID,
                            design=design,
                            sexID=sexID,
                            celltypeID=celltypeID,
                            y=y,
                            region=region,
                            control=control,
                            pval_adjust_method=pval_adjust_method,
                            rmv_zero_count_genes=rmv_zero_count_genes,
                            coeff=coeff,
                            fdr=fdr,
                            nom_pval=nom_pval,
                            Nperms=Nperms)
    # create power plots
    power_plots(SCE=SCE,
                range_downsampled=range_downsampled_cells,
                output_path=output_path,
                sampled="cells",
                sampleID=sampleID,
                celltypeID=celltypeID,
                fdr=fdr,
                nom_pval=nom_pval,
                Nperms=Nperms)
    # create down-sampling correlation plots (individuals)
    downsampling_corrplots(SCE=SCE,
                           range_downsampled=range_downsampled_individuals,
                           output_path=output_path,
                           sampled="individuals",
                           celltypeID=celltypeID,
                           Nperms=Nperms)
    # create down-sampling correlation plots (cells)
    downsampling_corrplots(SCE=SCE,
                           range_downsampled=range_downsampled_cells,
                           output_path=output_path,
                           sampled="cells",
                           celltypeID=celltypeID,
                           Nperms=Nperms)
}