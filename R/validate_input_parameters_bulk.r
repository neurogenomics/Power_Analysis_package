#' Tests input parameters for functions

#' @param SCEs A list of SingleCellExperiment (SCE) objects, each representing a scRNA-seq dataset.
#' @param dataset_names A vector of names corresponding to each dataset (as you would like them to appear in output plots).
#' @param celltype the cell type we are focusing on (name as it appears in cell type sub-directory name)
#' @param celltype_correspondence A named vector that maps a standard cell type label (e.g., `"Endo"`, `"Micro"`) to how that cell type appears in each dataset. Use `NA` if the cell type is not present in a given dataset.
#' @param output_path A clean directory path where down-sampled outputs and plots will be saved (should contain no subdirectories).
#' @param range_downsampled vector or list containing values which the data will be downsampled at, in ascending order
#' @param celltypeIDs A character vector specifying the column name in each SCE that denotes cell type identity (in order of SCEs).
#' @param sampleIDs  A character vector specifying the column name in each SCE that represents sample or donor IDs (in order of SCEs).
#' @param sampled Specifies the unit of down-sampling. Can be either `"individuals"` or `"cells"`, depending on whether the analysis downsamples across samples or cells.
#' @param bulkDE DGE analysis output for a bulk RNA-seq dataset (e.g., `LFSR.tsv`): rows (rownames) should be the genes, columns should be tissues, and entries should be significance levels
#' @param bulk_cutoff Proportion (0–1) of bulk tissues in which a gene must be differentially expressed to be considered (e.g., 0.9 selects DEGs found in ≥90% of tissues).
#' @param pvalue P-value threshold for defining DEGs in the bulk dataset.
#' @param Nperms Number of permutations to perform for each down-sampling level. Default is 20.
#' @param fontsize_axislabels Font size for axis labels in plot
#' @param fontsize_axisticks Font size for axis tick labels in plot
#' @param fontsize_title Font size for plot title
#' @param fontsize_legendlabels Font size for legend labels in plot
#' @param fontsize_legendtitle Font size for legend title in plot
#' @param plot_title Plot title

#' Checks all bulk analysis parameters are specified correctly

validate_input_parameters_bulk <- function(SCEs="placeholder",
                                           dataset_names="placeholder",
                                           celltype="placeholder",
                                           celltype_correspondence="placeholder",
                                           output_path="placeholder",
                                           range_downsampled="placeholder",
                                           celltypeIDs="placeholder",
                                           sampled="placeholder",
                                           sampleIDs="placeholder",
                                           bulkDE="placeholder",
                                           bulk_cutoff="placeholder",
                                           pvalue="placeholder",
                                           Nperms="placeholder",
                                           fontsize_axislabels="placeholder",
                                           fontsize_axisticks="placeholder",
                                           fontsize_title="placeholder",
                                           fontsize_legendlabels="placeholder",
                                           fontsize_legendtitle="placeholder",
                                           plot_title="placeholder"){

    # test each parameter to check if it works
    if(!identical(SCEs, "placeholder")){
        if(!is.list(SCEs)){
            stop("Error: SCEs should be a list of SingleCellExperiment objects.")
        }
    }
    if(!identical(dataset_names, "placeholder")){
        if(!is.character(dataset_names)&!is.list(dataset_names)){
            stop("Error: dataset_names should be a list of strings specifying the names of the datasets.")
        }
    }
    if(celltype!="placeholder"){
        if(!is.character(celltype)){
            stop("Error: celltype should be a string specifying the cell type.")
        }
    }
    if (!identical(celltype_correspondence, "placeholder")) {
        if (!is.list(celltype_correspondence)) {
            stop("Error: celltype_correspondence should be a list of lists containing the names of cell types as they appear across all DGE analysis outputs.")
        }
        for (celltype_list in celltype_correspondence) {
            if (!is.character(celltype_list)) {
                stop("Error: each element of celltype_correspondence should be a vector of cell type names.")
            }
            if(!identical(SCEs, "placeholder")){
                if(length(celltype_list)!=length(SCEs)){
                    stop("Error: each element of celltype_correspondence should have the same length as the number of datasets in SCEs.")
                }
            }
        }
    }
    if(output_path!="placeholder"){
        if(output_path!=getwd()){
            if(!is.character(output_path)){
                stop("Error: output_path should be a string specifying the base path where down-sampled DE analysis outputs are saved (and output will be saved).")
            }
        }
        if(!dir.exists(output_path)){
            stop("Error: the specified output_path directory does not exist.")
        }
        if(length(list.dirs(output_path, full.names = FALSE, recursive = FALSE)) > 1){
            stop(paste0("Error: The specified output_path directory is not clean (contains subdirectories): '", output_path, "'. Please provide an empty directory."))
        }
    }
    if(!identical(range_downsampled,"placeholder")){
        if(class(range_downsampled)=="character"){
            stop("Error: range_downsampled should be a list or numeric vector containing the values to be downsampled at.")
        }else{
            if(!is.list(range_downsampled)&!is.numeric(range_downsampled)){
                stop("Error: range_downsampled should be a list or numeric vector containing the values to be downsampled at.")
            }else if(is.unsorted(unlist(range_downsampled))|!sum(1*(range_downsampled > 0))==length(range_downsampled)){
                stop("Error: range_downsampled should be an ascending list with all elements bigger than 0.")
            }
        }
    }
    if(!identical(celltypeIDs,"placeholder")){
        if(!is.character(celltypeIDs)&!is.list(celltypeIDs)){
            stop("Error: celltypeIDs should be a string or list/vector specifying the cell type IDs in order of SCEs.")
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
    if(!identical(sampleIDs,"placeholder")){
        if(!is.character(sampleIDs)&!is.list(sampleIDs)){
            stop("Error: sampleIDs should be a string or list/vector specifying the cell type IDs in order of SCEs.")
        }
    }
    if(!identical(bulkDE,"placeholder")){
        if(class(bulkDE)!="data.frame"){
            stop("Error: bulkDE should be a dataframe with rows being genes, columns being tissues and entries being significance level")
        }
    }
    if(bulk_cutoff!="placeholder"){
        if(class(bulk_cutoff)!="numeric"){
            stop("Error: bulk_cutoff should be numerical.")
        }else{
            if(bulk_cutoff<0|bulk_cutoff>1){
                stop("Error: bulk_cutoff should be between 0 and 1.")
            }
        }
    }
    if(pvalue!="placeholder"){
        if(class(pvalue)!="numeric"){
            stop("Error: pvalue should be numerical.")
        }else{
            if(pvalue<0|pvalue>1){
                stop("Error: pvalue should be between 0 and 1.")
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
    if(fontsize_axisticks!="placeholder"){
        if(class(fontsize_axisticks)!="numeric"){
            stop("Error: fontsize_axisticks should be numerical.")
        }else{
            if(fontsize_axisticks-floor(fontsize_axisticks)!=0|fontsize_axisticks<0){
                stop("Error: fontsize_axisticks should be a positive integer.")
            }
        }
    }
    if(fontsize_title!="placeholder"){
        if(class(fontsize_title)!="numeric"){
            stop("Error: fontsize_title should be numerical.")
        }else{
            if(fontsize_title-floor(fontsize_title)!=0|fontsize_title<0){
                stop("Error: fontsize_title should be a positive integer.")
            }
        }
    }
    if(fontsize_legendlabels!="placeholder"){
        if(class(fontsize_legendlabels)!="numeric"){
            stop("Error: fontsize_legendlabels should be numerical.")
        }else{
            if(fontsize_legendlabels-floor(fontsize_legendlabels)!=0|fontsize_legendlabels<0){
                stop("Error: fontsize_legendlabels should be a positive integer.")
            }
        }
    }
    if(fontsize_legendtitle!="placeholder"){
        if(class(fontsize_legendtitle)!="numeric"){
            stop("Error: fontsize_legendtitle should be numerical.")
        }else{
            if(fontsize_legendtitle-floor(fontsize_legendtitle)!=0|fontsize_legendtitle<0){
                stop("Error: fontsize_legendtitle should be a positive integer.")
            }
        }
    }
    if(plot_title!="placeholder"){
        if(class(plot_title)!="character"){
            stop("Error: plot_title should be text.")
        }
    }

}
