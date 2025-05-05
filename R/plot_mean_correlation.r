#' Obtain the average correlation (across celltypes) at a specified cutoff p-value

#' @importFrom ggcorrplot ggcorrplot
#' @importFrom ggplot2 theme element_text
#' @importFrom utils write.csv

#' @param main_dataset name of the dataset used to select significant DEGs from (specified as a string, name as in dataset_names)
#' @param SCEs list of the input data (elements should be SCE objects)
#' @param sampleIDs list or vector of sample IDs (in order of SCEs)
#' @param celltypeIDs list or vector of cell type IDs (in order of SCEs)
#' @param pvals the cut-off p-value which will be used to select DEGs
#' @param celltype_correspondence list of different names specifying each cell type
#' @param dataset_names names of the datasets as they appear in the correlation plot (in order of SCEs)
#' @param sex_DEGs If TRUE, only keep genes present on sex chromosmomes. Queries hspanies gene Ensembl dataset.
#' @param output_path base path in which outputs will be stored

#' Saves mean correlation matrix (and actual values) in the appropriate directory

plot_mean_correlation <- function(main_dataset,
                                  SCEs,
                                  sampleIDs,
                                  celltypeIDs,
                                  pvals,
                                  celltype_correspondence,
                                  dataset_names,
                                  sex_DEGs=FALSE,
                                  output_path=getwd()){

    # validate function input params
    validate_input_parameters_correlation(main_dataset=main_dataset, SCEs=SCEs, sampleIDs=sampleIDs, celltypeIDs=celltypeIDs, pvalues=pvals,
                                          celltype_correspondence=celltype_correspondence, dataset_names=dataset_names, sex_DEGs=sex_DEGs,
                                          output_path=output_path)
    # outputs
    output_list <- list()

    # get DEouts
    DEouts <- list()
    for(idx in seq_along(SCEs)){
        dataset <- SCEs[[idx]]
        coeff_use <- as.character(sort(unique(colData(dataset)$sex))[[2]])
        savepath <- file.path(output_path,dataset_names[[idx]])
        DEouts[[idx]] <- DGE_analysis(SCE=dataset, design=~sex, sampleID=sampleIDs[[idx]], celltypeID=celltypeIDs[[idx]], coef=coeff_use, output_path=savepath)
    }

    # loop over each p-value
    for(pvalue in pvals){
        # list for genes of each celltype at specified p-value
        genes <- list()
        allCorrs <- list()
        i <- 0
        for(celltype in names(celltype_correspondence)){
            # get corresponding cell type names for each dataset
            celltype_names <- celltype_correspondence[[celltype]]
            # correlation for each celltype at specified p-value
            corrOut <- plot_celltype_correlation(main_dataset, DEouts, celltype_names, dataset_names, sex_DEGs, pvalue)
            i <- i+1
            # get correlation matrix for each celltype
            allCorrs[[i]] <- corrOut[[1]]
            # get all present genes for current celltype
            genes[[i]] <- list(corrOut[[4]])
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

        # plot correlation matrix
        corr_plot.plot <- ggcorrplot(round(meanCorr,3),
        hc.order = F,insig="pch",pch=5,pch.col = "grey",
        pch.cex=9,
        title=paste0("Total ", totNumGenes," DEGs"),
        colors = c("#FC4E07", "white", "#00AFBB"),
        outline.color = "white", lab = TRUE, lab_size=3.5,
        sig.level=0.05) + theme(plot.title = element_text(hjust = 0.7)) # add/remove type="upper" in ggcorrplot (after hc.order) to get upper triangular/full matrix

        # store output in list with p-value as key
        output_list[[as.character(pvalue)]] <- list(corr_plot = corr_plot.plot, meanCorr = meanCorr)

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
