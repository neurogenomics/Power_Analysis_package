#' Obtain highly randomised permutations of a specified dataset (based on sex labels), while maintaining consistency with the sample ID

#' @importFrom SingleCellExperiment colData
#' @importFrom infotheo mutinformation

#' @param SCE the input data (should be an SCE object)
#' @param sampleID sample ID
#' @param sexID sex ID
#' @param Nrandom_perms the number of randomised permutations the user needs (default 5)
#' @param Ntests the number of random permutation tests (default 100)

#' @return a list of size Nrandom_perms, with each item being a randomised dataset

random_permutations <- function(SCE,
                                sampleID="Donor.ID",
                                sexID="sex",
                                Nrandom_perms=5,
                                Ntests=100){

    # check input parameters are fine
	if(class(SCE)[1]!="SingleCellExperiment"){
        stop("Error: SCE should be a SingleCellExperiment object.")
	}
	if(sampleID!="Donor.ID"){	
		if(!is.character(sampleID)){
			stop("Error: sampleID should be a string specifying the column name.")
		}
		if(is.null(SCE[[sampleID]])){
			stop("Error: The specified column name (sampleID argument) is not present in the SCE. Please check if this is spelled correctly.")
		}
	}
	if(Nrandom_perms!=5){
		if(!as.integer(Nrandom_perms)==Nrandom_perms){
			stop("Error: Nrandom_perms should be a positive whole number.")
		}
	}
	if(Ntests!=100){	
		if(!as.integer(Ntests)==Ntests){
			stop("Error: Ntests should be a positive whole number.")
		}
	}	

    # output
    out <- list()

    # dataframe containing all of the samples and sexes
    uniq_key <- unique(colData(SCE)[,c(sampleID,sexID)])
    # dataframe containing all samples and sexes, at the individual cell level, so including all cells
    key <- colData(SCE)[,c(sampleID,sexID)]
    key$rand_sample_id <- key[[sampleID]]
    # list of subjects
    subjs <- uniq_key[[sampleID]]
    set.seed(101)

    for(i in 1:Ntests){
        # full data
        rand_data <- SCE
        # randomise sample IDs for each permutation
        uniq_key$rand_sample_id <- sample(uniq_key[[sampleID]])
        # set new sex labels for full dataset
        currentPerm <- merge(uniq_key,key[,c("rand_sample_id",sexID)], by="rand_sample_id")
        # select and rename relevant columns
        currentPerm <- currentPerm[c("rand_sample_id",paste0(sexID,".x"))]
        names(currentPerm) <- c("rand_sample_id","sex")
        # alter entire dataset
        rand_data[[sampleID]] <- currentPerm$rand_sample_id
        rand_data[[sexID]] <- currentPerm$sex
        # store in list
        out[[i]] <- rand_data
    }
    # pick top Nrandom_perms "most different to data"
    # sort in same order as perms
    data_new <- SCE
    orderedSubjs <- unique(colData(out[[1]])[[sampleID]])
    key_subjs <- data.frame(key=orderedSubjs,weight=1:length(orderedSubjs))
    merged <- merge(colData(SCE),key_subjs,by.x=sampleID,by.y='key',all.x=T,all.y=F)
    res <- merged[order(merged$weight),c(sampleID,sexID)]
    data_new[[sampleID]] <- res[[sampleID]]
    data_new[[sexID]] <- res[[sexID]]
    # now get MI and pick best Nrandom_perms
    orderedPerms <- list()
    for(i in 1:Ntests){
        orderedPerms[[i]] <- mutinformation(colData(out[[i]])[[sexID]], colData(data_new)[[sexID]])
    }
    # sort and get indices of best Nrandom_perms
    bestIndices <- order(unlist(orderedPerms))[1:Nrandom_perms]
    bestPerms <- list()
    for(j in 1:Nrandom_perms){
        bestPerms[[j]] <- out[[bestIndices[[j]]]]
    }

    return(bestPerms)

}