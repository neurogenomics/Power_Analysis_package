# Define global variables
utils::globalVariables(c(".","name"))

#' Computes the correlation of the log-foldchange of the DEGs (at specified p-value) for a given celltype, across all datasets

#' @importFrom ggcorrplot ggcorrplot
#' @importFrom data.table setkey rbindlist
#' @importFrom stats reshape complete.cases cor

#' @param main_dataset name of the dataset used to select significant DEGs from (specified as a string, name as in dataset_names)
#' @param DEouts a list containing outputs of DGE analysis (as returned/optionally saved by DGE_analysis) for datasets to be used in the correlation analysis
#' @param celltypes_list list of different names specifying each cell type (in order for each dataset in DEouts)
#' @param dataset_names names of the datasets as they will appear in the correlation plot
#' @param sex_DEGs If TRUE, only keep genes present on sex chromosmomes. Queries hspanies gene Ensembl dataset.
#' @param pval the cut-off p-value which will be used to select DEGs (default is 1 to include all the genes)

#' @return correlation matrix, plot and the number of DEGs at the specified p-value

plot_celltype_correlation <- function(main_dataset,
                                      DEouts,
                                      celltypes_list,
                                      dataset_names,
                                      sex_DEGs=FALSE,
                                      pval=1){

    # check input parameters are fine
    if(!is.character(main_dataset)){
        stop("Error: main_dataset should be a string specifying the dataset to select significant DEGs from")
    }
    if(class(DEouts)!="list"){
        stop("Error: DEouts should be a list.")
    }
    if(!is.list(celltypes_list) & !is.vector(celltypes_list)){
        stop("Error: celltypes_list should be a vector or list containing the names of the cell type of concern, as they appear across each DGE analysis output (in order).")
    }
    if(class(sex_DEGs)!="logical"){
        stop("Error: sex_DEGs should be TRUE (if DEGs are chosen only from sex chromosomes) or FALSE")
    }
    if(!is.numeric(pval) | pval <= 0 | pval > 1){
        stop("Error: pval should be a (positive) number between 0 and 1.")
    }
    if(!is.vector(dataset_names)){
        stop("Error: dataset_names should be a vector or list of strings specifying the names of the datasets as they should appear in the correlation plot.")
    }

    # filter sex_DEGs
    if (sex_DEGs) {
        DEouts <- lapply(DEouts, function(x) {
            x$celltype_all_genes <- sex_chromosome_DEGs(x$celltype_all_genes)
            # Only keep DEGs for genes present in filtered all_genes
            celltype_DEGs <- lapply(names(x$celltype_DEGs), function(y) {
                celltype_DEG <- x$celltype_DEGs[[y]][x$celltype_DEGs[[y]]$name %in% x$celltype_all_genes[[y]]$name,]
                return(celltype_DEG)
            })
            names(celltype_DEGs) <- names(x$celltype_DEGs)
            x$celltype_DEGs <- celltype_DEGs
            return(x)
        })
    }

    # variable to redefine DEouts
    allstudies <- list()
    # for each study, select data corresponding only to celltype
    for(i in seq_along(DEouts)){
        # get corresponding cell type names for current study
        celltype_name <- celltypes_list[[i]][[1]]
        # redefine DEouts so each element only contains study/celltype
        allstudies[[i]] <- DEouts[[i]]$celltype_all_genes[[celltype_name]]
    }
    # get names of datasets
    names(allstudies) <- dataset_names

    # reshape data so "dataset" is now a variable
    allstudies_dt <- rbindlist(allstudies,idcol="dataset")
    # create list of all genes across datasets
    allGenes <- list()
    for(j in 1:length(dataset_names)){
        allGenes[[j]] <- allstudies_dt[dataset==dataset_names[[j]],name]
    }
    # get set of genes in all studies
    shared_genes <- Reduce(intersect,allGenes)
    # number of [celltype] genes shared across all studies
    #print(length(shared_genes))
    #filter dataset (to only include shared genes)
    setkey(allstudies_dt,name)
    allstudies_dt <- allstudies_dt[shared_genes]

    # select DEGs by specifying only significant genes
    genes <- allstudies_dt[dataset==main_dataset & PValue<pval, name]
    # filter to just these genes
    allstudies_dt <- allstudies_dt[name %in% genes,]

    # make matrix for corr() - cols will be [datset, name, logFC]
    mat_lfc <- allstudies_dt[,.(dataset,name,logFC)]
    # reshape so cols now are [(gene) name, logFC.dataset1, logFC.dataset2,...,logFC.datasetN] (N being length(names(DEouts)))
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
    print(paste0("Number of ",celltype_name," genes at given p value is ",num_genes))
    genes <- mat_lfc$name # (only updates if we remove NA genes, else is the same as "genes" above)

    # get correlation matrix
    corr_lfc <- cor(mat_lfc[,2:ncol(mat_lfc)],method = "spearman") #### could also try "pearson"...

    # plot correlation matrix
    corr_plot.plot <- ggcorrplot(round(corr_lfc,2),
            hc.order = F, insig="pch",pch=5,pch.col = "grey",
            pch.cex=9,
            title=paste0("LFC Correlation Matrix, ",celltype_name," - pval < ",pval,": ",
                        num_genes," genes"),
            colors = c("#FC4E07", "white", "#00AFBB"),
            outline.color = "white", lab = TRUE,
            sig.level=0.05) # add/remove type="upper" in ggcorrplot (after hc.order) to get upper triangular/full matrix

    # output matrix, plot and genecount
    return(list(corr_lfc, corr_plot.plot, num_genes, genes))

}
