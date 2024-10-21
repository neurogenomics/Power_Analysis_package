# Define global variables
utils::globalVariables(c(".","name"))

#' Computes the correlation of the log-foldchange of the DEGs (at specified p-value) for a given celltype, across all datasets

#' @importFrom ggcorrplot ggcorrplot
#' @importFrom data.table setkey rbindlist
#' @importFrom stats reshape complete.cases cor

#' @param dataset_name name of the dataset used to select significant DEGs from (specified as a string, name as in allstudies)
#' @param allstudies a list containing all the datasets (most likely as SCE objects)
#' @param celltype the celltype to compute correlation for
#' @param p_val the cut-off p-value which will be used to select DEGs (default is 1 to include all the genes)

#' @return correlation matrix, plot and the number of DEGs at the specified p-value

plot_celltype_correlation <- function(dataset_name,
                                      allstudies,
                                      celltype,
                                      p_val=1){

    # check input parameters are fine
    if(!is.character(dataset_name)){
        stop("Error: dataset_name should be a string specifying the dataset to select significant DEGs from")
    }
    if(class(allstudies)!="list"){
        stop("Error: allstudies should be a list.")
    }
    if(!is.character(celltype)){
        stop("Error: celltype should be a string specifying the cell type to compute correlation for.")
    }
    if(!is.numeric(p_val) | p_val <= 0 | p_val > 1){
        stop("Error: p_val should be a (positive) number between 0 and 1.")
    }

    # variable to redefine allstudies
    allstudies_new <- list()
    j <- 0
    # for each study, select data corresponding only to celltype
    for(study in names(allstudies)){
        j <- j+1
        # redefine allstudies so each element only contains study/celltype
        allstudies_new[[j]] <- allstudies[[study]][[celltype]]
        names(allstudies_new)[[j]] <- study
    }

    # reshape data so "dataset" is now a variable
    allstudies_dt <- rbindlist(allstudies_new,idcol="dataset")
    # get names of datasets
    datasets <- unique(allstudies_dt$dataset)
    # create list of all genes across datasets
    allGenes <- list()
    for(j in 1:length(datasets)){
        allGenes[[j]] <- allstudies_dt[dataset==datasets[[j]],name]
    }
    # get set of genes in all studies
    shared_genes <- Reduce(intersect,allGenes)  
    # number of [celltype] genes shared across all studies
    #print(length(shared_genes))
    #filter dataset (to only include shared genes)
    setkey(allstudies_dt,name)
    allstudies_dt <- allstudies_dt[shared_genes]

    # select DEGs by specifying only significant genes
    genes <- allstudies_dt[dataset==dataset_name & PValue<p_val, name]
    # filter to just these genes
    allstudies_dt <- allstudies_dt[name %in% genes,]
    
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
    # remove NA rows // convert NA to 0s
    mat_lfc <- mat_lfc[complete.cases(mat_lfc), ]
    #### either above line (remove NA) OR below 2 (NAs to 0)
    #mat_lfc <- as.data.frame(mat_lfc)
    #mat_lfc[is.na(mat_lfc)] <- 0

    # get number of DEGs for this celltype
    num_genes <- dim(mat_lfc)[[1]]
    print(paste0("Number of ",celltype," genes at given p value is ",num_genes))
    genes <- mat_lfc$name # (only updates if we remove NA genes, else is the same as "genes" above)

    # get correlation matrix
    corr_lfc <- cor(mat_lfc[,2:ncol(mat_lfc)],method = "spearman") #### could also try "pearson"...

    # plot correlation matrix
    corr_plot.plot <- ggcorrplot(round(corr_lfc,2), 
            hc.order = F, insig="pch",pch=5,pch.col = "grey",
            pch.cex=9,
            title=paste0("LFC Correlation Matrix, ",celltype," - pval < ",p_val,": ",
                        num_genes," genes"),
            colors = c("#FC4E07", "white", "#00AFBB"),
            outline.color = "white", lab = TRUE,
            sig.level=0.05) # add/remove type="upper" in ggcorrplot (after hc.order) to get upper triangular/full matrix

    # output matrix, plot and genecount
    return(list(corr_lfc, corr_plot.plot, num_genes, genes))

}
