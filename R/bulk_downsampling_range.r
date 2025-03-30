#' Obtain the range of values to downsample at (either for individuals, or mean number of cells per individual), for bulk analysis

#' @importFrom SingleCellExperiment colData

#' @param SCEs list of the input data (elements should be SCE objects)
#' @param celltype_correspondence list of different names specifying each cell type
#' @param sampled downsampling carried out based on what (either "individuals" or "cells")
#' @param sampleIDs list or vector of sample IDs (in order of SCEs)
#' @param celltypeIDs list or vector of cell type IDs (in order of SCEs)

#' @return list containing values which the data will be downsampled at, in ascending order

bulk_downsampling_range <- function(SCEs,
                                    celltype_correspondence,
                                    sampled="individuals",
                                    sampleIDs="donor_id",
                                    celltypeIDs="cell_type"){
    
    # check sampleIDs, celltypeIDs
    if(length(sampleIDs) == 1){
        sampleIDs <- rep(sampleIDs,length(SCEs))
    }
    if(length(celltypeIDs) == 1){
        celltypeIDs <- rep(celltypeIDs,length(SCEs))
    }
    # initialise list to hold samples, largest dataset
    max_samples <- 0
    largest_dataset <- NULL
    # loop through cell type mapping
    for(standard_celltype in names(celltype_correspondence)){
        # loop through datasets to get downsampled range
        for(idx in seq_along(SCEs)){
            dataset <- SCEs[[idx]]
            celltype_name <- celltype_correspondence[[standard_celltype]][[idx]]
            if(!is.na(celltype_name)){
                # subset dataset
                dataset1 <- dataset[, colData(dataset)[[celltypeIDs[[idx]]]] == celltype_name]
                num_samples <- length(unique(colData(dataset1)[[sampleIDs[[idx]]]]))
                # check if this dataset/celltype has the most samples
                if(num_samples > max_samples){
                    max_samples <- num_samples
                    largest_dataset <- dataset1
                    sampleID <- sampleIDs[[idx]]
                }
            }
        }
    }
    
    # get downsampling range, return
    range_downsampled <- downsampling_range(largest_dataset,sampled,sampleID)
    return(range_downsampled)

}