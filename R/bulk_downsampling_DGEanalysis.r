#' Downsample the dataset, based either on the individuals or cells, and run DE analysis on each downsampled output. Save results in a dataframe

#' @importFrom SingleCellExperiment colData

#' @param SCEs list of the input data (elements should be SCE objects)
#' @param dataset_names list of the names of the datasets (as you would like them to appear in the "output_path" directory)
#' @param celltype_correspondence list of different names specifying each cell type
#' @param sampled downsampling carried out based on what (either "individuals" or "cells")
#' @param sampleIDs list or vector of sample IDs (in order of SCEs)
#' @param celltypeIDs list or vector of cell type IDs (in order of SCEs)
#' @param output_path base output directory where down-sampled DGE analysis outputs will be saved

#' Saves DGE analysis output in the correct directory, to be used by other bulk analysis functions

bulk_downsampling_DGEanalysis <- function(SCEs,
                                          dataset_names,
                                          celltype_correspondence,
                                          sampled="individuals",
                                          sampleIDs="donor_id",
                                          celltypeIDs="cell_type",
                                          output_path=getwd()){
                                    
    # validate function input params
    validate_input_parameters_bulk(SCEs=SCEs, dataset_names=dataset_names, celltype_correspondence=celltype_correspondence, sampled=sampled,
                                   sampleIDs=sampleIDs, celltypeIDs=celltypeIDs, output_path=output_path)

    # get biggest downsampling range
    max_downsampling_range <- bulk_downsampling_range(SCEs=SCEs, sampled=sampled, sampleIDs=sampleIDs, celltypeIDs=celltypeIDs, celltype_correspondence=celltype_correspondence)
    # check sampleIDs, celltypeIDs
    if(length(sampleIDs) == 1){
        sampleIDs <- rep(sampleIDs,length(SCEs))
    }
    if(length(celltypeIDs) == 1){
        celltypeIDs <- rep(celltypeIDs,length(SCEs))
    }
    # loop through all datasets and cell types
    for(standard_celltype in names(celltype_correspondence)){
        for(idx in seq_along(SCEs)){
            dataset <- SCEs[[idx]]
            celltype_name <- celltype_correspondence[[standard_celltype]][[idx]]
            if(!is.na(celltype_name)){
                # subset dataset                
                dataset1 <- dataset[, colData(dataset)[[celltypeIDs[[idx]]]] == celltype_name]
                # run downsampling, DGE analysis
                numsamples <- length(unique(colData(dataset1)[[sampleIDs[[idx]]]]))
                range_dataset <- max_downsampling_range[max_downsampling_range <= numsamples]
                savepath <- file.path(output_path,dataset_names[[idx]],standard_celltype)
                # run
                downsampling_DEanalysis(dataset1,range_dataset,output_path=savepath,sampleID=sampleIDs[[idx]],celltypeID=celltypeIDs[[idx]],coeff="M")
            }
        }
    }
}