#' Tests input parameters for functions

#' @param main_dataset name of the dataset used to select significant DEGs from (specified as a string, name as in dataset_names)
#' @param SCEs list of the input data (elements should be SCE objects)
#' @param sampleIDs list or vector of sample IDs (in order of SCEs)
#' @param celltypeIDs list or vector of cell type IDs (in order of SCEs)
#' @param celltype_correspondence list of different names specifying each cell type
#' @param pvalues the list of cut-off p-values which will be used to select DEGs (can just provide a list with one as well)
#' @param dataset_names names of the datasets as they appear in the correlation plot (in order of SCEs)
#' @param corr_mats (named) list of correlation matrices for each celltype with the final element being the mean correlation matrix, all at specified p-value
#' @param num_real_datasets total number of *real* datasets (most likely the number of studies, but sometimes a study may be split e.g. into 2 brain regions, so in this case it would be the number of studies plus 1)
#' @param alphaval (alpha) transparency of the non-mean boxplots
#' @param N_randperms number of random permutations of the dataset used to select significant DEGs from
#' @param N_subsets number of pairs of random subsets of the dataset used to select significant DEGs from
#' @param sex_DEGs true if DEGs come from sex chromosomes, else false
#' @param fontsize_yaxislabels font size for axis labels in plot
#' @param fontsize_yaxisticks font size for axis tick labels in plot
#' @param fontsize_title font size for plot title
#' @param fontsize_legendlabels font size for legend labels in plot
#' @param fontsize_legendtitle font size for legend title in plot
#' @param fontsize_facet_labels font size for facet labels
#' @param output_path base path in which outputs will be stored

#' Checks all correlation analysis parameters are specified correctly

validate_input_parameters_correlation <- function(main_dataset="placeholder",
                                                  SCEs="placeholder",
                                                  sampleIDs="placeholder",
                                                  celltypeIDs="placeholder",
                                                  celltype_correspondence="placeholder",
                                                  pvalues="placeholder",
                                                  dataset_names="placeholder",
                                                  corr_mats="placeholder",
                                                  num_real_datasets="placeholder",
                                                  alphaval="placeholder",
                                                  N_randperms="placeholder",
                                                  N_subsets="placeholder",
                                                  sex_DEGs="placeholder",
                                                  fontsize_yaxislabels="placeholder",
                                                  fontsize_yaxisticks="placeholder",
                                                  fontsize_title="placeholder",
                                                  fontsize_legendlabels="placeholder",
                                                  fontsize_legendtitle="placeholder",
                                                  fontsize_facet_labels="placeholder",
                                                  output_path="placeholder"){

    # test each parameter to check if it works
    if(main_dataset!="placeholder"){
        if(!is.character(main_dataset)){
            stop("Error: main_dataset should be a string")
        }
    }
    if(!identical(SCEs, "placeholder")){
        if(!is.list(SCEs)){
            stop("Error: SCEs should be a list of SingleCellExperiment objects.")
        }
    }
    if(!identical(sampleIDs,"placeholder")){
        if(!is.character(sampleIDs)&!is.list(sampleIDs)){
            stop("Error: sampleIDs should be a string or list/vector specifying the cell type IDs in order of SCEs.")
        }
    }
    if(!identical(celltypeIDs,"placeholder")){
        if(!is.character(celltypeIDs)&!is.list(celltypeIDs)){
            stop("Error: celltypeIDs should be a string or list/vector specifying the cell type IDs in order of SCEs.")
        }
    }
    if(!identical(celltype_correspondence,"placeholder")){    
        if(class(celltype_correspondence)!="list"){
            stop("Error: celltype_correspondence should be a list of lists of all cell type names, as they appear in each of the SCEs")
        }
    }
    if(!identical(pvalues,"placeholder")){
        for(pval in pvalues){
            if(!is.numeric(pval) | pval <= 0 | pval > 1){
                stop("Error: pvalues should be a list of (positive) numbers between 0 and 1")
            }
        }
    }
    if(!identical(dataset_names,"placeholder") && !is.vector(dataset_names)){
        stop("Error: dataset_names should be a list containing the names of all datasets as should appear in the final output (if these are different to the names in SCEs)")
    }
    if(!identical(corr_mats,"placeholder")){
        if(class(corr_mats)!="list"){
            stop("Error: corr_mats should be a list of matrices")
        }
    }
    if(num_real_datasets!="placeholder"){
        if(floor(num_real_datasets)!=num_real_datasets){
            stop("Error: num_real_datasets should be an integer specifying the number of real datasets (see description)")
        }
    }
    if(alphaval!="placeholder"){
        if(class(alphaval)!="numeric"){
            stop("Error: alphaval should be numerical")
        }
    }
    if(N_randperms!="placeholder"){
        if(floor(N_randperms)!=N_randperms){
            stop("Error: N_randperms should be an integer specifying the number of random permutations of Tsai (see description)")
        }
    }
    if(N_subsets!="placeholder"){
        if(floor(N_subsets)!=N_subsets){
            stop("Error: N_subsets should be an integer specifying the number of subsets of Tsai (see description)")
        }
    }
    if(sex_DEGs!="placeholder"){
        if(class(sex_DEGs)!="logical"){
            stop("Error: sex_DEGs should be TRUE (if DEGs are chosen only from sex chromosomes) or FALSE")
        }
    }
    if(fontsize_yaxislabels!="placeholder"){
        if(class(fontsize_yaxislabels)!="numeric"){
            stop("Error: fontsize_yaxislabels should be numerical")
        }else{
            if(fontsize_yaxislabels-floor(fontsize_yaxislabels)!=0|fontsize_yaxislabels<0){
                stop("Error: fontsize_yaxislabels should be a positive integer")
            }
        }
    }
    if(fontsize_yaxisticks!="placeholder"){
        if(class(fontsize_yaxisticks)!="numeric"){
            stop("Error: fontsize_yaxisticks should be numerical")
        }else{
            if(fontsize_yaxisticks-floor(fontsize_yaxisticks)!=0|fontsize_yaxisticks<0){
                stop("Error: fontsize_yaxisticks should be a positive integer")
            }
        }
    }
    if(fontsize_title!="placeholder"){
        if(class(fontsize_title)!="numeric"){
            stop("Error: fontsize_title should be numerical")
        }else{
            if(fontsize_title-floor(fontsize_title)!=0|fontsize_title<0){
                stop("Error: fontsize_title should be a positive integer")
            }
        }
    }
    if(fontsize_legendlabels!="placeholder"){
        if(class(fontsize_legendlabels)!="numeric"){
            stop("Error: fontsize_legendlabels should be numerical")
        }else{
            if(fontsize_legendlabels-floor(fontsize_legendlabels)!=0|fontsize_legendlabels<0){
                stop("Error: fontsize_legendlabels should be a positive integer")
            }
        }
    }
    if(fontsize_legendtitle!="placeholder"){
        if(class(fontsize_legendtitle)!="numeric"){
            stop("Error: fontsize_legendtitle should be numerical")
        }else{
            if(fontsize_legendtitle-floor(fontsize_legendtitle)!=0|fontsize_legendtitle<0){
                stop("Error: fontsize_legendtitle should be a positive integer")
            }
        }
    }
    if(fontsize_facet_labels!="placeholder"){
        if(class(fontsize_facet_labels)!="numeric"){
            stop("Error: fontsize_facet_labels should be numerical")
        }else{
            if(fontsize_facet_labels-floor(fontsize_facet_labels)!=0|fontsize_facet_labels<0){
                stop("Error: fontsize_facet_labels should be a positive integer")
            }
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
}