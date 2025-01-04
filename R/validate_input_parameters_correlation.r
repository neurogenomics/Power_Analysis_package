#' Tests input parameters for functions

#' @param dataset_name name of the dataset used to select significant DEGs from (specified as a string, name as in allStudies)
#' @param allstudies a list containing all the datasets (most likely as SCE objects)
#' @param celltypes a list containing the celltypes to compute mean correlation across
#' @param pvalue the cut-off p-value which will be used to select DEGs
#' @param data_names names of the datasets as they appear in the correlation plot
#' @param corrMats (named) list of correlation matrices for each celltype with the final element being the mean correlation matrix, all at specified p-value
#' @param numRealDatasets total number of *real* datasets (most likely the number of studies, but sometimes a study may be split e.g. into 2 brain regions, so in this case it would be the number of studies plus 1)
#' @param alphaval (alpha) transparency of the non-mean boxplots
#' @param numPerms number of random permutations of the dataset used to select significant DEGs from
#' @param numSubsets number of pairs of random subsets of the dataset used to select significant DEGs from
#' @param sexDEGs true if DEGs come from sex chromosomes, else false
#' @param fontsize_yaxislabels font size for axis labels in plot
#' @param fontsize_yaxisticks font size for axis tick labels in plot
#' @param fontsize_title font size for plot title
#' @param fontsize_legendlabels font size for legend labels in plot
#' @param fontsize_legendtitle font size for legend title in plot
#' @param fontsize_facet_labels font size for facet labels
#' @param output_path base path in which outputs will be stored

#' Checks all correlation analysis parameters are specified correctly

validate_input_parameters_correlation <- function(dataset_name="placeholder",
                                                  allstudies="placeholder",
                                                  celltypes="placeholder",
                                                  pvalue="placeholder",
                                                  data_names="placeholder",
                                                  corrMats="placeholder",
                                                  numRealDatasets="placeholder",
                                                  alphaval="placeholder",
                                                  numPerms="placeholder",
                                                  numSubsets="placeholder",
                                                  sexDEGs="placeholder",
                                                  fontsize_yaxislabels="placeholder",
                                                  fontsize_yaxisticks="placeholder",
                                                  fontsize_title="placeholder",
                                                  fontsize_legendlabels="placeholder",
                                                  fontsize_legendtitle="placeholder",
                                                  fontsize_facet_labels="placeholder",
                                                  output_path="placeholder"){

    # test each parameter to check if it works
    if(dataset_name!="placeholder"){
        if(!is.character(dataset_name)){
            stop("Error: dataset_name should be a string")
        }
    }
    if(allstudies!="placeholder"){
        if(class(allstudies)!="list"){
            stop("Error: allstudies should be a list")
        }
    }
    if(celltypes!="placeholder"){    
        if(!is.character(celltypes)){
            stop("Error: celltypes should be a list containing strings specifying the celltypes")
        }
    }
    if(pvalue!="placeholder"){
        if(!is.numeric(pvalue) | pvalue <= 0 | pvalue > 1){
            stop("Error: pvalue should be a (positive) number between 0 and 1")
        }
    }
    if(data_names!="placeholder" && !is.vector(data_names)){
        stop("Error: data_names should be a list containing the names of all datasets as should appear in the final output (if these are different to the names in allstudies)")
    }
    if(corrMats!="placeholder"){
        if(class(corrMats)!="list"){
            stop("Error: corrMats should be a list of matrices")
        }
    }
    if(numRealDatasets!="placeholder"){
        if(floor(numRealDatasets)!=numRealDatasets){
            stop("Error: numRealDatasets should be an integer specifying the number of real datasets (see description)")
        }
    }
    if(alphaval!="placeholder"){
        if(class(alphaval)!="numeric"){
            stop("Error: alphaval should be numerical")
        }
    }
    if(numPerms!="placeholder"){
        if(floor(numPerms)!=numPerms){
            stop("Error: numPerms should be an integer specifying the number of random permutations of Tsai (see description)")
        }
    }
    if(numSubsets!="placeholder"){
        if(floor(numSubsets)!=numSubsets){
            stop("Error: numSubsets should be an integer specifying the number of subsets of Tsai (see description)")
        }
    }
    if(sexDEGs!="placeholder"){
        if(class(sexDEGs)!="logical"){
            stop("Error: sexDEGs should be TRUE (if DEGs are chosen only from sex chromosomes) or FALSE")
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