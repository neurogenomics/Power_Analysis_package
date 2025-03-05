#' Obtain the range of values to downsample at (either for individuals, or mean number of cells per individual), for bulk analysis

#' @importFrom SingleCellExperiment colData

#' @param datasets list of the input data (elements should be SCE objects)
#' @param celltype_correspondence list of different names specifying each cell type
#' @param sampled downsampling carried out based on what (either "individuals" or "cells")
#' @param sampleID sample ID
#' @param celltype_ID cell type ID

#' @return list containing values which the data will be downsampled at, in ascending order

bulk_downsampling_range <- function(datasets,
                                    celltype_correspondence,
                                    sampled="individuals",
                                    sampleID="donor_id",
                                    celltype_ID="cell_type"){

    # initialise list to hold samples, largest dataset
    max_samples <- 0
    largest_dataset <- NULL
    # loop through cell type mapping
    for(standard_celltype in names(celltype_correspondence)){
        # loop through datasets to get downsampled range
        for(idx in seq_along(datasets)){
            dataset <- datasets[[idx]]
            celltype_name <- celltype_correspondence[[standard_celltype]][[idx]]
            dataset1 <- dataset[, colData(dataset)[[celltype_ID]] == celltype_name]
            num_samples <- length(unique(colData(dataset1)[[sampleID]]))
            # check if this dataset/celltype has the most samples
            if(num_samples > max_samples){
                max_samples <- num_samples
                largest_dataset <- dataset1
            }
        }
    }
    # get downsampling range, return
    range_downsampled <- downsampling_range(largest_dataset,sampled,sampleID)
    return(range_downsampled)

}