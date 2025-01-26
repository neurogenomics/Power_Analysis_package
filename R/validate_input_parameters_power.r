#' Tests input parameters for functions

#' @param data the input data (should be an SCE object)
#' @param range_downsampled vector or list containing values which the data will be downsampled at, in ascending order
#' @param output_path base path in which outputs will be stored
#' @param inpath base path where downsampled DGE analysis output is stored (taken to be output_path if not provided)
#' @param sampled downsampling carried out based on what (either "individuals" or "cells")
#' @param sampleID sample ID
#' @param design the design formula of class type `formula`. Equation used to fit the model- data for the generalised linear model e.g. expression ~ sex + pmi + disease
#' @param sexID sex ID
#' @param celltypeID cell type ID
#' @param coeff which coefficient to carry out DE analysis with respect to
#' @param fdr the cut-off False Discovery Rate below which to select DEGs
#' @param nom_pval the cut-off nominal P-value below which to select DEGs (as an alternative to FDR)
#' @param Nperms number of subsets created when downsampling at each level
#' @param N_randperms number of randomised permutations of the dataset (based on sex) to be correlated
#' @param N_subsetpairs number of pairs of subsets of the dataset to be correlated
#' @param y the column name in the SCE object for the return variable e.g. "diagnosis" - Case or disease. Default is the last variable in the design formula. y can be discrete (logistic regression) or continuous (linear regression)
#' @param region the column name in the SCE object for the study region. If there are multiple regions in the study (for example two brain regions). Pseudobulk values can be derived separately. Default is "single_region" which will not split by region.
#' @param control character specifying which control level for the differential expression analysis e.g. in a case/control/other study use "control" in the y column to compare against. NOTE only need to specify if more than two groups in y, leave as default value for two groups or continuous y. Default is NULL.
#' @param pval_adjust_method the adjustment method for the p-value in the differential expression analysis. Default is benjamini hochberg "BH". See  stats::p.adjust for available options
#' @param rmv_zero_count_genes whether genes with no count values in any cell should be removed. Default is TRUE
#' @param base_outpath the base output path which contains power analysis output folders for *all* datasets (each one as created by "downsampling_DEanalysis.r")
#' @param datasets_path path where all dataset SCEs (in .qs format) are saved
#' @param sampleIDs list of sample ID names for datasets in "datasets_path" (may be in any order)

#' Checks all power analysis parameters are specified correctly

validate_input_parameters_power <- function(data="placeholder",
                                            range_downsampled="placeholder",
                                            output_path="placeholder",
                                            inpath="placeholder",
                                            sampled="placeholder",
                                            sampleID="placeholder",
                                            design="placeholder",
                                            sexID="placeholder",
                                            celltypeID="placeholder",
                                            coeff="placeholder",
                                            fdr="placeholder",
                                            nom_pval="placeholder",
                                            Nperms="placeholder",
                                            N_randperms="placeholder",
                                            N_subsetpairs="placeholder",
                                            y="placeholder",
                                            region="placeholder",
                                            control="placeholder",
                                            pval_adjust_method="placeholder",
                                            rmv_zero_count_genes="placeholder",
                                            base_outpath="placeholder",
                                            datasets_path="placeholder",
                                            sampleIDs="placeholder"){

    # test each parameter to check if it works
    if(class(data)[1]!="SingleCellExperiment"){
        stop("Error: data should be a SingleCellExperiment object.")
    }
    if(class(range_downsampled)=="character"){
        if(!identical(range_downsampled,"placeholder")){
            stop("Error: range_downsampled should be a list or numeric vector containing the values to be downsampled at.")
        }
    }else{
        if(!is.list(range_downsampled)&!is.numeric(range_downsampled)){
            stop("Error: range_downsampled should be a list or numeric vector containing the values to be downsampled at.")
        }else if(is.unsorted(unlist(range_downsampled))|!sum(1*(range_downsampled > 0))==length(range_downsampled)){
            stop("Error: range_downsampled should be an ascending list with all elements bigger than 0.")
        }else if(sum(unlist(range_downsampled)-floor(unlist(range_downsampled)))!=0){
            stop("Error: range_downsampled should contain only integer values.")
        }
    }
    if(output_path!="placeholder"){
        if(output_path!=getwd()){
            if(!is.character(output_path)){
                stop("Error: output_path should be a string specifying the base path where output will be stored.")        
            }
            if(!dir.exists(output_path)){
                stop("Error: the specified output_path directory does not exist.")
            }      
        }
    }
    if(inpath!="placeholder"){
        if(inpath!=output_path){
            if(!is.character(inpath)){
                stop("Error: inpath should be a string specifying the base path where down-sampling DE analysis outputs are saved.")        
            }
            if(!dir.exists(inpath)){
                stop("Error: the specified inpath directory does not exist.")
            }      
        }
    }
    if(sampled!="placeholder"){
        if(sampled!="individuals"){
            if(!is.character(sampled)){
                stop("Error: sampled should be a string specifying what we are downsampling based on.")
            }
            if(sampled!="cells"){
                stop("Error: sampled should either be set to individuals or cells.")
            }
        }
    }
    if(sampleID!="placeholder"){
        if(sampleID!="Donor.ID"){	
            if(!is.character(sampleID)){
                stop("Error: sampleID should be a string specifying the column name.")
            }
        }
    }
    if(design!="placeholder"){
        if(design!=as.formula(paste0("~",sexID))){
            if(class(design)!="formula"){
                stop("Error: design should be a formula (e.g. ~ var1+var2+...)")
            }
        }
    }
    if(sexID!="placeholder"){
        if(sexID!="sex"){	
            if(!is.character(sexID)){
                stop("Error: sexID should be a string specifying the column name.")
            }
            if(is.null(data[[sexID]])){
                stop("Error: The specified column name (sexID argument) is not present in the data. Please check if this is spelled correctly.")
            }
        }
    }
    if(celltypeID!="placeholder"){
        if(celltypeID!="cell_type"){	
            if(!is.character(celltypeID)){
                stop("Error: celltypeID should be a string specifying the column name.")
            }
            if(is.null(data[[celltypeID]])){
                stop("Error: The specified column name (celltypeID argument) is not present in the data. Please check if this is spelled correctly.")
            }
        }
    }
    if(coeff!="placeholder"){
        if(coeff!="male"){	
            if(!is.character(coeff)){
                stop("Error: coeff should be a string specifying coefficient DE analysis should be carried out with respect to.")
            }
        }
    }
    if(fdr!="placeholder"){
        if(class(fdr)!="numeric"){
            stop("Error: fdr should be numerical.")
        }else{
            if(fdr<0|fdr>1){
                stop("Error: fdr should be between 0 and 1.")
            }
        }
    }
    if(nom_pval!="placeholder"){
        if(class(nom_pval)!="numeric"){
            stop("Error: nom_pval should be numerical.")
        }else{
            if(nom_pval<0|nom_pval>1){
                stop("Error: nom_pval should be between 0 and 1.")
            }
        }
    }
    if(Nperms!="placeholder"){
        if(class(Nperms)!="numeric"){
            stop("Error: Nperms should be numerical.")
        }else{
            if(Nperms < 0|Nperms%%1!=0){
                stop("Error: Nperms should be a positive integer.")
            }
        }
    }
    if(N_randperms!="placeholder"){
        if(class(N_randperms)!="numeric"){
            stop("Error: N_randperms should be numerical.")
        }else{
            if(N_randperms < 0|N_randperms%%1!=0){
                stop("Error: N_randperms should be a positive integer.")
            }
        }
    }
    if(N_subsetpairs!="placeholder"){
        if(class(N_subsetpairs)!="numeric"){
            stop("Error: N_subsetpairs should be numerical.")
        }else{
            if(N_subsetpairs < 0|N_subsetpairs%%1!=0){
                stop("Error: N_subsetpairs should be a positive integer.")
            }
        }
    }
    if(!is.null(y)){
        if(y!="placeholder"){
            if(!is.character(y)){
                stop("Error: please input a character for y indicating the column holding the response variable information")
            }
            if(is.null(data[[y]])){
                stop("Error: the inputted y value is not present in the data, perhaps check the spelling is correct")
            }
        }
    }
    if(region!="placeholder"){
        if(region!="single_region"){
            if(!is.character(region)){
                stop("Error: please input a character for region indicating the column holding the variable information")
            }
            if(is.null(data[[region]])){
                stop("Error: the inputted region value is not present in the data, perhaps check the spelling is correct")
            }
        }
    }
    if(!is.null(control)){
        if(control!="placeholder"){
            if(!control %in% unique(data[[y]])){
                stop("Error: the inputted control value is not present in y in the data, perhaps check the spelling is correct")
            }
        }
    }
    if(pval_adjust_method!="placeholder"){
        if(!is.character(pval_adjust_method)){
            stop("Error: please input a character for the pval_adjust_method indicating the method to be used")
        }
    }
    if(rmv_zero_count_genes!="placeholder"){
        if(!is.logical(rmv_zero_count_genes)){
            stop("Error: please input TRUE/FALSE for rmv_zero_count_genes")
        }
    }
    if(base_outpath!="placeholder"){
        if(length(list.dirs(base_outpath,recursive=F))!=0){
            if(!is.character(base_outpath)){
                stop("Error: base_outpath should be a string specifying the base path where all power analysis outputs are stored.")        
            }
            if(!dir.exists(base_outpath)){
                stop("Error: the specified base_outpath directory does not exist.")
            }      
        }else{
            stop("Error: base_outpath cannot be an empty directory.")
        }
    }
    if(datasets_path!="placeholder"){
        if(length(list.files(datasets_path))!=0){
            if(!is.character(datasets_path)){
                stop("Error: datasets_path should be a string specifying the path where all datasets are stored.")        
            }
            if(!dir.exists(datasets_path)){
                stop("Error: the specified datasets_path directory does not exist.")
            }      
        }else{
            stop("Error: datasets_path cannot be an empty directory.")
        }
    }
    if(class(sampleIDs)!="character"&!is.list(sampleIDs)){
        stop("Error: sampleIDs should be a list or vector containing the list of sampleID names for each dataset.")
    }
}