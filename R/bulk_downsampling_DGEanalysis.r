#' Downsample the dataset, based either on the individuals or cells, and run DE analysis on each downsampled output. Save results in a dataframe

#' @importFrom SingleCellExperiment colData

#' @param SCEs A list of SingleCellExperiment (SCE) objects, each representing a scRNA-seq dataset.
#' @param dataset_names A vector of names corresponding to each dataset (as you would like them to appear in output plots).
#' @param celltype_correspondence A named vector that maps a standard cell type label (e.g., `"Endo"`, `"Micro"`) to how that cell type appears in each dataset. Use `NA` if the cell type is not present in a given dataset.
#' @param output_path A directory path where down-sampled outputs and plots will be saved.
#' @param celltypeIDs A character vector specifying the column name in each SCE that denotes cell type identity (in order of SCEs).
#' @param sampleIDs  A character vector specifying the column name in each SCE that represents sample or donor IDs (in order of SCEs).
#' @param sampled Specifies the unit of down-sampling. Can be either `"individuals"` or `"cells"`, depending on whether the analysis downsamples across samples or cells.
#' @param pvalue the cut-off p-value used to select DEGs
#' @param Nperms number of permutations of DGE analysis outputs for each sample

#' Saves DGE analysis output in the correct directory, to be used by other bulk analysis functions

bulk_downsampling_DGEanalysis <- function(SCEs,
                                          dataset_names,
                                          celltype_correspondence,
                                          sampled = c("individuals", "cells"),
                                          sampleIDs="donor_id",
                                          celltypeIDs="cell_type",
                                          output_path=getwd(),
                                          pvalue=0.05,
                                          Nperms=20){

    # validate function input params
    validate_input_parameters_bulk(SCEs=SCEs, dataset_names=dataset_names, celltype_correspondence=celltype_correspondence, sampled=sampled,
                                   sampleIDs=sampleIDs, celltypeIDs=celltypeIDs, output_path=output_path, pvalue=pvalue, Nperms=Nperms)

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
            coeff_use <- as.character(sort(unique(colData(dataset)$sex))[[2]])
            celltype_name <- celltype_correspondence[[standard_celltype]][[idx]]
            if(!is.na(celltype_name)){
                # subset dataset
                dataset1 <- dataset[, colData(dataset)[[celltypeIDs[[idx]]]] == celltype_name]
                # run downsampling, DGE analysis
                numsamples <- length(unique(colData(dataset1)[[sampleIDs[[idx]]]]))
                range_dataset <- max_downsampling_range[max_downsampling_range <= numsamples]
                savepath <- file.path(output_path,dataset_names[[idx]],standard_celltype)
                # run
                downsampling_DEanalysis(dataset1,range_dataset,output_path=savepath,sampled=sampled,
                                        sampleID=sampleIDs[[idx]],celltypeID=celltypeIDs[[idx]],coeff=coeff_use,
                                        nom_pval=pvalue,Nperms=Nperms)
            }
        }
    }
}
