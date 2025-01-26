#' Obtain the range of values to downsample at (either for individuals, or mean number of cells per individual)

#' @importFrom SingleCellExperiment colData
#' @importFrom stats quantile

#' @param data the input data (should be an SCE object)
#' @param sampled downsampling carried out based on what (either "individuals" or "cells")
#' @param sampleID sample ID

#' @return list containing values which the data will be downsampled at, in ascending order

downsampling_range <- function(data,
                               sampled="individuals",
                               sampleID="Donor.ID"){

    # get list of sample IDs
    coldata <- colData(data)
    IDs <- unique(coldata[[sampleID]])
    # get number of samples
    numSamples <- length(IDs)

    if(sampled=="individuals"){
        
        # want to pick 10 (equally spaced out) sampling points
        range_downsampled <- round(ceiling(seq(from=0,to=numSamples,length.out=12))/5)*5
        range_downsampled <- unique(range_downsampled[2:11])

    }else{

        # get number of cells per individual
        numCells <- list()
        for(i in 1:numSamples){
            numCells[[i]] <- dim(coldata[coldata[[sampleID]]==IDs[[i]],])[[1]]
        }
        numCells <- unlist(numCells)
        # pick 10 sampling points
        range_downsampled <- round(seq(0,unname(round(quantile(numCells,0.9)/5)*5), length.out=11)/5)*5
        range_downsampled <- unique(range_downsampled[2:11])

    }

    return(range_downsampled)

} 