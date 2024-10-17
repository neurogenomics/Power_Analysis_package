# Define global variables
utils::globalVariables(c(".","dataset"))

#' For a given down-sampled DE output, computes the correlation of the log-foldchange of the DEGs (at specified p-value) for a given dataset (celltype)

#' @importFrom ggcorrplot ggcorrplot
#' @importFrom data.table rbindlist setkey
#' @importFrom stats reshape complete.cases cor
#' @importFrom ggplot2 theme element_text scale_fill_gradient2 labs

#' @param allstudies - a list containing DGE analysis outputs (subset to "celltype_all_genes" and the specified cell type) at each down-sampled point
#' @param sampled - downsampling carried out based on what (either "individuals" or "cells")

#' @return correlation matrix, plot and the number of DEGs at the specified p-value

compute_downsampled_corr <- function(allstudies,
                                     sampled="individuals"){

    # check input parameters are fine
    if(class(allstudies)!="list"){
        stop("Error: allstudies should be a list.")
    }

    # variable to redefine allstudies
    allstudies_new <- list()
    j <- 1
    # select data for each down-sampling point
    for(sampling_point in names(allstudies)){
        # redefine allstudies so each element only contains data for that down-sampled point
        allstudies_new[[j]] <- allstudies[[sampling_point]]
        names(allstudies_new)[[j]] <- sampling_point
        j <- j+1
    }
    # reshape data so "dataset" is now a variable
    allstudies_dt <- rbindlist(allstudies_new,idcol="dataset")
    #filter dataset
    setkey(allstudies_dt,name)
    
    # make matrix for corr() - cols will be [datset, name, logFC]
    mat_lfc <- allstudies_dt[,.(dataset,name,logFC)]
    # reshape so cols now are [(gene) name, logFC.dataset1, logFC.dataset2,...,logFC.datasetN] (N being length(names(allstudies)))
    mat_lfc <-
    reshape(mat_lfc, idvar = "name", timevar = "dataset", 
            direction = "wide")
    # remove logFC. from name (now cols are just [name, dataset1, dataset2,...,datasetN] with logFC values in each column)
    colnames(mat_lfc)<-
        c(colnames(mat_lfc)[1],
            substr(colnames(mat_lfc[,2:ncol(mat_lfc)]),
                    7,nchar(colnames(mat_lfc[,2:ncol(mat_lfc)]))))
    # remove NA rows
    mat_lfc <- mat_lfc[complete.cases(mat_lfc), ]

    # get correlation matrix
    corr_lfc <- cor(mat_lfc[,2:ncol(mat_lfc)],method = "spearman")

    # plot correlation matrix
    if(sampled=="individuals"){
        corr_plot.plot <- ggcorrplot(round(corr_lfc,2), 
                hc.order = F, insig="pch",pch=5,pch.col = "grey",
                pch.cex=9,
                title=paste0("LogFC Correlation Matrix when downsampling individuals"),
                colors = c("#FC4E07", "white", "#00AFBB"),
                outline.color = "white", lab = TRUE, lab_size=3.5,
                sig.level=0.05) + 
                theme(plot.title = element_text(hjust = 0.7)) + 
                scale_fill_gradient2(breaks=c(0, 1), limit=c(0, 1)) +
                labs(fill="Corr")
    }else{
        corr_plot.plot <- ggcorrplot(round(corr_lfc,2), 
        hc.order = F, insig="pch",pch=5,pch.col = "grey",
        pch.cex=9,
        title=paste0("LogFC Correlation Matrix when downsampling cells"),
        colors = c("#FC4E07", "white", "#00AFBB"),
        outline.color = "white", lab = TRUE, lab_size=3.5,
        sig.level=0.05) + 
        theme(plot.title = element_text(hjust = 0.7)) +
        scale_fill_gradient2(breaks=c(0, 1), limit=c(0, 1)) +
        labs(fill="Corr")
    }

    # output matrix, plot
    return(list(corr_lfc, corr_plot.plot))

}