#' Validate that there are no user input errors for the differential expression analysis - sc.cell.type.de

#' @param SCE SingleCellExperiment object, a specialised S4 class for storing data from single-cell experiments
#' @param design Design formula of class type `formula`. Equation used to fit the model- data for the generalised linear model.
#' @param sampleID Column name in the SCE object to perform pseudobulk on, usually the patient identifier. This column is used for grouping in the pseudobulk approach
#' @param celltypeID Column name in the SCE object for the cell type variable 
#' @param y Column name in the SCE object for the return variable e.g. "diagnosis" - Case or disease
#' @param region Column name in the SCE object for the study region. If there are multiple regions in the study (for example two brain regions). Pseudobulk values can be derived separately. Default is "single_region" which will not split by region.
#' @param coef Character specifying which level to investigate for the differential expression analysis e.g. in a case/control study use "case" if case is the identifier in the y column to get positive fold changes to relate to case samples. leave as default value for continuous y. Default is NULL.
#' @param control Character specifying which control level for the differential expression analysis e.g. in a case/control/other study use "control" in the y column to compare against. NOTE only need to specify if more than two groups in y, leave as default value for two groups or continuous y. Default is NULL.
#' @param pval_adjust_method Adjustment method for the p-value in the differential expression analysis. Default is benjamini hochberg "BH". See  stats::p.adjust for available options
#' @param adj_pval Adjusted p-value cut-off for the differential expression analysis, 0-1 range
#' @param output_path Folder where the graphs from the differential expression analysis are saved. Default will create a folder in the current working directory "sc.cell.type.de.graphs". False will skip plotting.
#' @param rmv_zero_count_genes Whether genes with no count values in any cell should be removed. Default is TRUE
#' @param verbose Logical indicating if extra information about the differential expression analysis should be printed

#' Checks all DE analysis parameters are specified correctly

validate_input_parameters_de<-function(SCE, design, sampleID, celltypeID, 
                                       y, region, coef, control, 
                                       pval_adjust_method, adj_pval, output_path, 
                                       rmv_zero_count_genes, verbose){
    if(class(SCE)[1]!="SingleCellExperiment")
        stop("Error: please input the data as a SingleCellExperiment object")
    if(!is.character(celltypeID))
        stop("Error: please input a character for celltypeID indicating the column holding the cell type information")
    if(is.null(SCE[[celltypeID]]))
        stop(paste0("Error: the inputted celltypeID: ",celltypeID,
                    " is not present in the SCE object, perhaps check if spelling is correct"))
    # check design is a formula, formatting is correct and variables exist
    if(!inherits(design,"formula"))
        stop("Error: please input a formula for the design variable specifying the comparison. See examples.")
    design_txt <- paste0(deparse(design,width.cutoff = 500),collapse=',')
    # check if formula contains `~` to ensure user inputting correct format
    if(!grepl( "~", design_txt, fixed = TRUE))
        stop("Error: please input a correctly formated formula for the design variable containing a `~`. See examples.")
    design_txt <- gsub(".*~","",design_txt)
    design_txt <- gsub("^\\s+|\\s+$","",strsplit(design_txt, "[+]")[[1]])
    # check for duplicate entries
    if(length(design_txt)!=length(unique(design_txt)))
        stop("Error: there are duplicate entries in the design formula")
    # check each variable to see if they are in SCE object
    for(i in design_txt){
        if(is.null(SCE[[i]]))
            stop(paste0("Error: the inputted value: ",i,
                        " in the design formula is not present in the SCE object, perhaps check for spelling is correct"))
    }
    # only check y if user is setting it themselves
    if(!is.null(y)){
        if(!is.character(y))
            stop("Error: please input a character for y indicating the column holding the response variable information")
        if(is.null(SCE[[y]]))
            stop(paste0("Error: the inputted y value: ",y,
                        " is not present in the SCE object, perhaps check for spelling is correct"))
    }
    else{
        # if y not specified take last value in design matrix
        y <- design_txt[[length(design_txt)]]
        if(is.null(SCE[[y]]))
            stop(paste0("Error: the inputted y value taken from your formula (the last variable): ",y,
                        " is not present in the SCE object, perhaps check for spelling is correct"))
    }
    # only check region if user is setting it themselves
    if(region!="single_region"){
        if(!is.character(region))
            stop("Error: please input a character for region indicating the column holding the variable information")
        if(is.null(SCE[[region]]))
            stop(paste0("Error: the inputted region value: ",region,
                        " is not present in the SCE object, perhaps check for spelling is correct"))
    }
    if(!is.null(coef)){
        if(!coef %in% unique(SCE[[y]]))
            stop(paste0("Error: the inputted coef value: ",coef,
                        " is not present in y ", y,
                        " in the SCE object, perhaps check for spelling is correct"))
    }
    if(!is.null(control)){
        if(!control %in% unique(SCE[[y]]))
            stop(paste0("Error: the inputted control value: ",control,
                        " is not present in y ", y,
                        " in the SCE object, perhaps check for spelling is correct"))
    }
    if(!is.character(sampleID))
        stop("Error: please input a character for sampleID indicating the column holding the patient identifier")
    if(is.null(SCE[[sampleID]]))
        stop(paste0("Error: the inputted patient_ID value: ",sampleID,
                    " is not present in the SCE object, perhaps check for spelling is correct"))
    if(!is.character(pval_adjust_method))
        stop("Error: please input a character for the pval_adjust_method indicating the method to be used")
    if(class(adj_pval)[1]!="numeric")
        stop("Error: please input a 0-1 range number for the adj_pval indicating the cut off of significance for the differential expression analysis")
    if(adj_pval<0|adj_pval>1)
        stop("Error: please input a 0-1 range number for the adj_pval indicating the cut off of significance for the differential expression analysis")
    # if the user doesn't want to save plots from DE analysis they can pass in a false for the folder input parameter
    if(!isFALSE(output_path))
        if(!is.character(output_path))
            stop("Error: for the output_path input variable, please pass a directory for the plots to be saved in, go with the default or pass in FALSE stop plotting")
    if(!is.logical(verbose))
        stop("Error: please input TRUE/FALSE for verbose")
    if(!is.logical(rmv_zero_count_genes))
        stop("Error: please input TRUE/FALSE for rmv_zero_count_genes")
}