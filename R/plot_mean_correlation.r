#' Obtain the average correlation (across celltypes) at a specified cutoff p-value

#' @importFrom ggcorrplot ggcorrplot
#' @importFrom ggplot2 theme element_text
#' @importFrom utils write.csv

#' @param dataset_name name of the dataset used to select significant DEGs from (specified as a string, name as in DEouts)
#' @param DEouts a list containing outputs of DGE analysis (as returned/optionally saved by DGE_analysis) for datasets to be used in the correlation analysis
#' @param pvals the cut-off p-value which will be used to select DEGs
#' @param celltype_correspondence list of different names specifying each cell type
#' @param data_names names of the datasets as they appear in the correlation plot
#' @param output_path base path in which outputs will be stored

#' Saves mean correlation matrix (and actual values) in the appropriate directory

plot_mean_correlation <- function(dataset_name,
                                  DEouts,
                                  pvals,
                                  celltype_correspondence,
                                  data_names,
                                  output_path=getwd()){

    # validate function input params
    validate_input_parameters_correlation(dataset_name=dataset_name, DEouts=DEouts, pvalues=pvals,
                                          celltype_correspondence=celltype_correspondence, data_names=data_names, output_path=output_path)
    # outputs
    output_list <- list()

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
            corrOut <- plot_celltype_correlation(dataset_name, DEouts, celltype_names, data_names, pvalue)
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
        if(is.vector(data_names)){
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

}
