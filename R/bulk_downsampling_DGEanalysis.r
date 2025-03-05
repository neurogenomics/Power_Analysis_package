#' Downsample the dataset, based either on the individuals or cells, and run DE analysis on each downsampled output. Save results in a dataframe

#' @importFrom SingleCellExperiment colData

#' @param datasets list of the input data (elements should be SCE objects)
#' @param celltype_correspondence list of different names specifying each cell type
#' @param sampled downsampling carried out based on what (either "individuals" or "cells")
#' @param sampleID sample ID
#' @param celltype_ID cell type ID
#' @param path base output directory where down-sampled DGE analysis outputs will be saved

#' Saves DGE analysis output in the correct directory, to be used by other bulk analysis functions

bulk_downsampling_DGEanalysis <- function(datasets,
                                          celltype_correspondence,
                                          sampled="individuals",
                                          sampleID="donor_id",
                                          celltype_ID="cell_type",
                                          path="placeholder"){
                                    
    # validate function input params
    validate_input_parameters_bulk(datasets=datasets, celltype_correspondence=celltype_correspondence, sampled=sampled,
                                   sampleID=sampleID, celltype_ID=celltype_ID, path=path)

    # get biggest downsampling range
    max_downsampling_range <- bulk_downsampling_range(datasets=datasets, sampled=sampled, sampleID=sampleID, celltype_ID=celltype_ID, celltype_correspondence=celltype_correspondence)
    # loop through all datasets and cell types
    for(standard_celltype in names(celltype_correspondence)){
        for(idx in seq_along(datasets)){
            dataset <- datasets[[idx]]
            celltype_name <- celltype_correspondence[[standard_celltype]][[idx]]
            dataset1 <- dataset[, colData(dataset)[[celltype_ID]] == celltype_name]
            # run downsampling, DGE analysis
            numsamples <- length(unique(colData(dataset1)[[sampleID]]))
            range_dataset <- max_downsampling_range[max_downsampling_range <= numsamples]
            savepath <- paste0(path,"/",datasets[idx],"/",standard_celltype,"/")
            # run
            downsampling_DEanalysis(dataset1,range_dataset,output_path=savepath,sampleID=sampleID,celltypeID=celltype_ID,coeff="M")
        }

    }

}