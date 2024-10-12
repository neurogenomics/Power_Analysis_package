#' Sample the data for a specified mean number of cells (per sample) from the dataset

#' @importFrom SingleCellExperiment colData

#' @param data the input data (should be an SCE object)
#' @param meanCells the mean number of cells to be sampled from each individual
#' @param Nperms the number of subsets the user needs (default 20)
#' @param sampleID sample ID

#' @return a list of size Nperms, with each item being a subset of the input SCE containing the specified mean number of cells per sample

sample_cells <- function(data,
                         meanCells,
                         sampleID="Donor.ID",
                         Nperms=20){

    # check input parameters are fine
	if(class(data)[1]!="SingleCellExperiment"){
        stop("Error: data should be a SingleCellExperiment object.")
	}
    if(!as.integer(meanCells)==meanCells){
        stop("Error: meanCells should be a positive whole number.")
    }
	if(sampleID!="Donor.ID"){	
		if(!is.character(sampleID)){
			stop("Error: sampleID should be a string specifying the column name.")
		}
		if(is.null(data[[sampleID]])){
			stop("Error: The specified column name (sampleID argument) is not present in the data. Please check if this is spelled correctly.")
		}
	}
	if(Nperms!=20){
		if(!as.integer(Nperms)==Nperms){
			stop("Error: Nperms should be a positive whole number.")
		}
	}
    
    # store outputs
    newData <- list()
    # column data
    coldata <- colData(data)
    # check no. cells per individual
    samples <- unique(coldata[[sampleID]])
    # from each individual, want to sample an average of meanCells cells
    # if the total number of cells for an individual < meanCells, just sample all cells in that case
    # list of cell IDs
    set.seed(101)
    for(j in 1:Nperms){
        cells <- list()
        for(sample in samples){
            # list of cells for this sample
            cells_tmp <- rownames(coldata[coldata[,sampleID]==sample,])
            # check if total number of cells for this sample < meanCells
            err <- sample(-5:5,1)
            if(length(cells_tmp) < meanCells + err){
                # if so, add all cells to list
                cells <- unlist(c(cells,cells_tmp))
            }else{
                # else, add meanCells cells to list
                sampleCells <- sample(cells_tmp, meanCells + err)
                cells <- unlist(c(cells, sampleCells))
            }
        }
        # get updated permutations and add to list
        newData[[j]] <- data[, unique(colnames(data))%in%cells]
    }

    return(newData)

}