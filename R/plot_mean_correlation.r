#' Obtain the average correlation (across celltypes) at a specified cutoff p-value

#' @importFrom ggcorrplot ggcorrplot
#' @importFrom ggplot2 theme element_text
#' @importFrom utils write.csv

#' @param main_dataset name of the dataset used to select significant DEGs from, and for which random permutations and subset pairs will be created (to be specified as a string, name as in dataset_names)
#' @param SCEs list of the input data (elements should be SCE objects)
#' @param sampleIDs list or vector of sample IDs (in order of SCEs)
#' @param celltypeIDs list or vector of cell type IDs (in order of SCEs)
#' @param celltype_correspondence list of different names specifying each cell type
#' @param dataset_names names of the datasets as they appear in the correlation plot (in order of SCEs)
#' @param N_randperms number of random permutations of the dataset to be created
#' @param N_subsets number of pairs of random subsets of the dataset to be created
#' @param pvals the cut-off p-value which will be used to select DEGs
#' @param sex_DEGs If TRUE, only keep genes present on sex chromosmomes. Queries hspanies gene Ensembl dataset.
#' @param output_path base path in which outputs will be stored

#' Saves mean correlation matrix (and actual values) in the appropriate directory

plot_mean_correlation <- function(main_dataset,
                                  SCEs,
                                  sampleIDs,
                                  celltypeIDs,
                                  celltype_correspondence,
                                  dataset_names,
                                  N_randperms=5,
                                  N_subsets=5,
                                  pvals=c(0.05,0.025,0.01,0.001,0.0001),
                                  sex_DEGs=FALSE,
                                  output_path=getwd()){

    # validate function input params
    validate_input_parameters_correlation(main_dataset=main_dataset, SCEs=SCEs, sampleIDs=sampleIDs, celltypeIDs=celltypeIDs,
                                          celltype_correspondence=celltype_correspondence, dataset_names=dataset_names, N_randperms=N_randperms,
                                          N_subsets=N_subsets, pvalues=pvals, sex_DEGs=sex_DEGs, output_path=output_path)
    # outputs
    output_list <- list()

    # get DEouts
    DEouts <- list()    
    for(idx in seq_along(SCEs)){
        dataset <- SCEs[[idx]]
        coeff_use <- as.character(sort(unique(colData(dataset)$sex))[[2]])
        savepath <- file.path(output_path,dataset_names[[idx]])
        print(paste0("Running DGE analysis for ", dataset_names[[idx]], "..."))
        DEouts[[idx]] <- DGE_analysis(SCE=dataset, design=~sex, sampleID=sampleIDs[[idx]], celltypeID=celltypeIDs[[idx]], coef=coeff_use, output_path=savepath)
    }

    # get DE outputs for random permutations and subsets
    idx_main <- which(dataset_names==main_dataset)
    if(length(idx_main) == 0){
        stop("Error: main_dataset not found in dataset_names.")
    }
    # get random permutations
    if(N_randperms > 0){
        print(paste0("Creating ", N_randperms, " random permutations of the main dataset..."))
        rand_perms <- random_permutations(SCEs[[idx_main]], sampleID=sampleIDs[[idx_main]], Nrandom_perms=N_randperms)
        for(i in seq_along(rand_perms)){
            DEouts[[length(DEouts)+1]] <- DGE_analysis(SCE=rand_perms[[i]], design=~sex, sampleID=sampleIDs[[idx_main]], celltypeID=celltypeIDs[[idx_main]], coef=coeff_use, output_path=file.path(output_path, paste0(main_dataset,"_randperm_", i)))
        }
        dataset_names <- c(dataset_names, paste0(main_dataset, "_RandPerm_", seq_len(N_randperms)))
    }
    # get random subsets
    if(N_subsets > 0){
        print(paste0("Creating ", N_subsets, " independent subsets of the main dataset..."))
        rand_subsets <- subset_pairs(SCEs[[idx_main]], sampleID=sampleIDs[[idx_main]], Noutputs=N_subsets)
        for(i in seq_along(rand_subsets)){
            subset_1 <- rand_subsets[[i]][[1]]
            subset_2 <- rand_subsets[[i]][[2]]
            # run DGE analysis for each subset pair
            DEouts[[length(DEouts)+1]] <- DGE_analysis(SCE=subset_1, design=~sex, sampleID=sampleIDs[[idx_main]], celltypeID=celltypeIDs[[idx_main]], coef=coeff_use, output_path=file.path(output_path, paste0(main_dataset,"_randsubset_", i, "a")))
            DEouts[[length(DEouts)+1]] <- DGE_analysis(SCE=subset_2, design=~sex, sampleID=sampleIDs[[idx_main]], celltypeID=celltypeIDs[[idx_main]], coef=coeff_use, output_path=file.path(output_path, paste0(main_dataset,"_randsubset_", i, "b")))
            dataset_names <- c(dataset_names, paste0(main_dataset, "_RandSubset_", i, "a"), paste0(main_dataset, "_RandSubset_", i, "b"))
        }
    }

    # loop over each p-value
    for(pvalue in pvals){
        # list for genes of each celltype at specified p-value
        genes <- list()
        allCorrs <- list()
        i <- 1
        for(celltype in names(celltype_correspondence)){
            # get corresponding cell type names for each dataset
            celltype_names <- celltype_correspondence[[celltype]]
            # add random permutations and subsets to celltype names
            if(N_randperms > 0){
                celltype_names <- c(celltype_names, rep(celltype_names[[idx_main]], N_randperms))
            }
            if(N_subsets > 0){
                celltype_names <- c(celltype_names, rep(celltype_names[[idx_main]], N_subsets))
            }
            # correlation for each celltype at specified p-value
            corrOut <- plot_celltype_correlation(main_dataset, DEouts, celltype_names, dataset_names, sex_DEGs, pvalue)
            # get correlation matrix for each celltype
            allCorrs[[celltype]] <- corrOut[[1]]
            # get all present genes for current celltype
            genes[[i]] <- list(corrOut[[4]])
            i <- i+1
        }
        # total number of unique genes across celltypes
        genes <- unlist(genes)
        totNumGenes <- length(unique(genes))
        # average correlations actoss celltypes, for specified p-value
        meanCorr <- Reduce("+",allCorrs)/length(allCorrs)

        # rename columns and rows
        if(is.vector(dataset_names)){
            rownames(meanCorr) <- colnames(meanCorr) <- dataset_names
        }
        # add meanCorr to the list of all correlations
        allCorrs[["Mean"]] <- meanCorr

        # plot correlation matrix
        corr_plot.plot <- ggcorrplot(round(meanCorr,3),
        hc.order = F,insig="pch",pch=5,pch.col = "grey",
        pch.cex=9,
        title=paste0("Total ", totNumGenes," DEGs"),
        colors = c("#FC4E07", "white", "#00AFBB"),
        outline.color = "white", lab = TRUE, lab_size=3.5,
        sig.level=0.05) + theme(plot.title = element_text(hjust = 0.7)) # add/remove type="upper" in ggcorrplot (after hc.order) to get upper triangular/full matrix

        # store output in list with p-value as key
        output_list[[as.character(pvalue)]] <- allCorrs

        # save the mean correlation matrix
        write.csv(meanCorr, file.path(output_path, paste0("mean_correlation_matrix_p", pvalue, ".csv")))
        # save the plot
        if(length(DEouts) <= 5){
            fig_width <- 10
            fig_height <- 10
        }else{
            fig_width <- 2*length(DEouts)
            fig_height <- 2*length(DEouts)
        }
        ggsave(file.path(output_path, paste0("mean_correlation_p", pvalue, ".png")), corr_plot.plot, width=fig_width, height=fig_height, units="cm", bg="white")
        ggsave(file.path(output_path, paste0("mean_correlation_p", pvalue, ".pdf")), corr_plot.plot, width=fig_width, height=fig_height, units="cm", bg="white")
    }
    return(output_list)
}
