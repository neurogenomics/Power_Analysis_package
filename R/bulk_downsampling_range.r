#' Obtain the range of values to downsample at (either for individuals, or mean number of cells per individual), for bulk analysis

#' @importFrom SingleCellExperiment colData

#' @param SCEs A list of SingleCellExperiment (SCE) objects, each representing a scRNA-seq dataset.
#' @param celltype_correspondence A named vector that maps a standard cell type label (e.g., `"Endo"`, `"Micro"`) to how that cell type appears in each dataset. Use `NA` if the cell type is not present in a given dataset.
#' @param celltypeIDs A character vector specifying the column name in each SCE that denotes cell type identity (in order of SCEs).
#' @param sampleIDs  A character vector specifying the column name in each SCE that represents sample or donor IDs (in order of SCEs).
#' @param sampled Specifies the unit of down-sampling. Can be either `"individuals"` or `"cells"`, depending on whether the analysis downsamples across samples or cells.


#' @return list containing values which the data will be downsampled at, in ascending order

bulk_downsampling_range <- function(SCEs,
                                    celltype_correspondence,
                                    sampled="individuals",
                                    sampleIDs="donor_id",
                                    celltypeIDs="cell_type"){

    sampled <- match.arg(sampled, choices = c("individuals", "cells"))

    # check sampleIDs, celltypeIDs
    if(length(sampleIDs) == 1){
        sampleIDs <- rep(sampleIDs,length(SCEs))
    }
    if(length(celltypeIDs) == 1){
        celltypeIDs <- rep(celltypeIDs,length(SCEs))
    }

    # initialise list to hold maximum number of samples/cells, largest dataset
    max_number <- 0
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

                # get number of samples or cells of the sub-dataset
                if (sampled == "individuals") {
                    num_units <- length(unique(colData(dataset1)[[sampleIDs[[idx]]]]))
                } else {  # sampled == "cells"
                    num_units <- ncol(dataset1)
                }

                # check if this dataset/celltype has the most samples/cells
                if (num_units > max_number) {
                    max_number  <- num_units
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
