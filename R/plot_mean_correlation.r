#' Obtain the average correlation (across celltypes) at a specified cutoff p-value

#' @importFrom ggcorrplot ggcorrplot
#' @importFrom ggplot2 theme element_text

#' @param dataset_name name of the dataset used to select significant DEGs from (specified as a string, name as in allStudies)
#' @param allstudies a list containing all the datasets (most likely as SCE objects)
#' @param celltypes a list containing the celltypes to compute mean correlation across
#' @param pvalue the cut-off p-value which will be used to select DEGs
#' @param data_names names of the datasets as they appear in the correlation plot

#' @return mean correlation matrix

plot_mean_correlation <- function(dataset_name,
                                  allstudies,
                                  celltypes,
                                  pvalue,
                                  data_names="placeholder"){

    # validate function input params
    validate_input_parameters_correlation(dataset_name=dataset_name, allstudies=allstudies, celltypes=celltypes,
                                          pvalue=pvalue, data_names=data_names)

    # list for genes of each celltype at specified p-value
    genes <- list()
    allCorrs <- list()
    i <- 0
    for(celltype in celltypes){
        # correlation for each celltype at specified p-value
        corrOut <- plot_celltype_correlation(dataset_name, allstudies, celltype, pvalue)
        i <- i+1
        # get correlation matrix for each celltype
        allCorrs[[i]] <- corrOut[[1]]
        # get all present genes for current celltype
        genes[[i]] <- corrOut[[4]]
    }
    # total number of unique genes across celltypes
    genes <- unlist(genes)
    totNumGenes <- length(unique(genes))
    # average correlations actoss celltypes, for specified p-value
    meanCorr <- Reduce("+",allCorrs)/length(allCorrs)

    # rename columns and rows
    if(data_names!="placeholder"&&is.vector(data_names)){
        rownames(meanCorr) <- colnames(meanCorr) <- data_names
    }

    # plot correlation matrix
    corr_plot.plot <- ggcorrplot(round(meanCorr,3), 
    hc.order = F,insig="pch",pch=5,pch.col = "grey",
    pch.cex=9,
    title=paste0("Total ", totNumGenes," DEGs"),
    colors = c("#FC4E07", "white", "#00AFBB"),
    outline.color = "white", lab = TRUE, lab_size=3.5,
    sig.level=0.05) + theme(plot.title = element_text(hjust = 0.7)) # add/remove type="upper" in ggcorrplot (after hc.order) to get upper triangular/full matrix

    # output plot and final matrix
    return(list(corr_plot.plot,meanCorr))
}
