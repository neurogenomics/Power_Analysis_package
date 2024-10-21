#' Obtain independent pairs of subsets of a specified dataset, based on sample ID

#' @importFrom SingleCellExperiment colData

#' @param data the input data (should be an SCE object)
#' @param sampleID sample ID
#' @param Noutputs the number of randomised permutations the user needs (default 5)

#' @return a list of size Noutputs, with each item being a pair of independent datasets

subset_pairs <- function(data,
                         sampleID="Donor.ID", 
                         Noutputs=5){

    # check input parameters are fine
	if(class(data)[1]!="SingleCellExperiment"){
        stop("Error: data should be a SingleCellExperiment object.")
	}
	if(sampleID!="Donor.ID"){	
		if(!is.character(sampleID)){
			stop("Error: sampleID should be a string specifying the column name.")
		}
		if(is.null(data[[sampleID]])){
			stop("Error: The specified column name (sampleID argument) is not present in the data. Please check if this is spelled correctly.")
		}
	}
	if(Noutputs!=5){
		if(!as.integer(Noutputs)==Noutputs){
			stop("Error: Noutputs should be a positive whole number.")
		}
	}

    # store outputs
    newData <- list()
    # loop through
    set.seed(101)
    for(j in 1:Noutputs){
        print(paste0("Creating subset ",j))
        # list of samples
        allSamples <- unique(colData(data)[[sampleID]])
        # pick both halves of the samples and subset to these
        subSamples1 <- sample(allSamples, length(allSamples)*0.5)
        subSamples2 <- allSamples[!allSamples %in% subSamples1]
        # get updated permutations and add to list
        newData1 <- data[, data[[sampleID]]%in%subSamples1]
        newData2 <- data[, data[[sampleID]]%in%subSamples2]
        newData[[j]] <- list(newData1,newData2)
    }

    return(newData)

}