#' Sample the data for a specified number of individuals from the dataset

#' @importFrom SingleCellExperiment colData

#' @param data the input data (should be an SCE object)
#' @param Nsamples the number of samples to be included in each subset
#' @param sampleID the column name in the SCE object which is the patient identifier (usually just "sample_id")
#' @param sexID the column name in the SCE object which is the sex (usually just "sex" or "Sex")
#' @param Nperms the number of subsets the user needs (default 20)

#' @return a list of size Nperms, with each item being a subset of the input SCE containing the specified number of individuals

sample_individuals <- function(data,
                               Nsamples, 
                               sampleID="Donor.ID", 
                               sexID="sex",
                               Nperms=20){

    # check input parameters are fine
	if(class(data)[1]!="SingleCellExperiment"){
        stop("Error: data should be a SingleCellExperiment object.")
	}
    if(!as.integer(Nsamples)==Nsamples){
        stop("Error: Nsamples should be a positive whole number.")
    }
	if(sampleID!="Donor.ID"){	
		if(!is.character(sampleID)){
			stop("Error: sampleID should be a string specifying the column name.")
		}
		if(is.null(data[[sampleID]])){
			stop("Error: The specified column name (sampleID argument) is not present in the data. Please check if this is spelled correctly.")
		}
	}
    if(sexID!="sex"){	
        if(!is.character(sexID)){
            stop("Error: sexID should be a string specifying the column name.")
        }
        if(is.null(data[[sexID]])){
            stop("Error: The specified column name (sexID argument) is not present in the data. Please check if this is spelled correctly.")
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
    # list of samples
    allSamples <- unique(coldata[[sampleID]])  
    # loop through
    set.seed(101)
    for(j in 1:Nperms){
        # number of males and females
        numMales <- 0
        numFemales <- 0
        print(paste0("Create subset ",j))
        while(numMales==0|numFemales==0){
            # pick the number of required samples and subset to these
            subSamples <- sample(allSamples, Nsamples)
            # update number of males, females
            for(sample in subSamples){
                if(unique(coldata[coldata[,sampleID]==sample,][[sexID]]) == unique(data[[sexID]])[[2]]){
                    numMales <- numMales + 1
                }else{
                    numFemales <- numFemales + 1
                }
            }
        }
        # get updated permutations and add to list
        newData[[j]] <- data[, data[[sampleID]]%in%subSamples]
    }

    return(newData)

}