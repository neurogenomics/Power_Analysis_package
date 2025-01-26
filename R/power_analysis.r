# Define global variables
utils::globalVariables(c("dataset"))

#' Runs entire power analysis pipeline

#' @importFrom stats as.formula

#' @param data the input data (should be an SCE object)
#' @param range_downsampled_individuals vector or list containing values which the data will be downsampled at (for individuals), in ascending order
#' @param range_downsampled_cells vector or list containing values which the data will be downsampled at (for cells), in ascending order
#' @param output_path base path in which outputs will be stored
#' @param sampleID sample ID
#' @param design the design formula of class type `formula`. Equation used to fit the model- data for the generalised linear model e.g. expression ~ sex + pmi + disease
#' @param sexID sex ID
#' @param celltypeID cell type ID
#' @param coeff which coefficient to carry out DE analysis with respect to
#' @param fdr the cut-off False Discovery Rate below which to select DEGs
#' @param nom_pval the cut-off nominal P-value below which to select DEGs (as an alternative to FDR)
#' @param Nperms number of subsets created when downsampling at each level
#' @param y the column name in the SCE object for the return variable e.g. "diagnosis" - Case or disease. Default is the last variable in the design formula. y can be discrete (logistic regression) or continuous (linear regression)
#' @param region the column name in the SCE object for the study region. If there are multiple regions in the study (for example two brain regions). Pseudobulk values can be derived separately. Default is "single_region" which will not split by region.
#' @param control character specifying which control level for the differential expression analysis e.g. in a case/control/other study use "control" in the y column to compare against. NOTE only need to specify if more than two groups in y, leave as default value for two groups or continuous y. Default is NULL.
#' @param pval_adjust_method the adjustment method for the p-value in the differential expression analysis. Default is benjamini hochberg "BH". See  stats::p.adjust for available options
#' @param rmv_zero_count_genes whether genes with no count values in any cell should be removed. Default is TRUE

#' Saves all plots and DGE analysis outputs in the appropriate directories
#' @export

power_analysis <- function(data,
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
    
    setwd(output_path)

    # alter range_downsampled_individuals
    if(identical(range_downsampled_individuals,"placeholder")){
        range_downsampled_individuals <- downsampling_range(data, "individuals", sampleID)
    }
    # alter range_downsampled_cells
    if(identical(range_downsampled_cells,"placeholder")){
        range_downsampled_cells <- downsampling_range(data, "cells", sampleID)
    }
    # alter design
    if(design=="placeholder"){
        design=as.formula(paste0("~",sexID))
    }    

    # create preliminary plots
    preliminary_plots(data=data, 
                      output_path=output_path, 
                      sampleID=sampleID, 
                      design=design, 
                      sexID=sexID, 
                      celltypeID=celltypeID, 
                      coeff=coeff, 
                      fdr=fdr)
    
    ## create power plots for down-sampling individuals, cells
    # down-sample individuals and run DE analysis
    downsampling_DEanalysis(data=data,
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
    power_plots(data=data, 
                range_downsampled=range_downsampled_individuals, 
                output_path=output_path,
                sampled="individuals", 
                sampleID=sampleID, 
                celltypeID=celltypeID, 
                fdr=fdr, 
                nom_pval=nom_pval, 
                Nperms=Nperms)
    # down-sample cells and run DE analysis
    downsampling_DEanalysis(data=data,
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
    power_plots(data=data, 
                range_downsampled=range_downsampled_cells, 
                output_path=output_path,
                sampled="cells", 
                sampleID=sampleID, 
                celltypeID=celltypeID, 
                fdr=fdr, 
                nom_pval=nom_pval, 
                Nperms=Nperms)
    # create down-sampling correlation plots (individuals)
    downsampling_corrplots(data=data, 
                           range_downsampled=range_downsampled_individuals, 
                           output_path=output_path, 
                           sampled="individuals", 
                           celltypeID=celltypeID, 
                           Nperms=Nperms)
    # create down-sampling correlation plots (cells)
    downsampling_corrplots(data=data, 
                           range_downsampled=range_downsampled_cells, 
                           output_path=output_path, 
                           sampled="cells", 
                           celltypeID=celltypeID, 
                           Nperms=Nperms)
}