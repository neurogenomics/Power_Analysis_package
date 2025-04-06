# Define global variables
utils::globalVariables(c(".","chromosome","x.chromosome_name"))

#' Given a DGE analysis output, this outputs a subset of the analysis with the genes laying only on the sex chromosomes (X/Y)

#' @importFrom data.table rbindlist setkey setorder :=
#' @importFrom biomaRt useDataset getBM useMart

#' @param all_genes a list whose elements are DGE analysis output for each celltype in a dataset (e.g. as returned by "DGE_analysis.R", but subsetted to "celltype_all_genes")
#' @param ensemblID_colname the column name within each DGE analysis output in all_genes, which contains the ensembl IDs for each gene

#' @return the DGE analysis output restricted to just the genes that lay on the sex chromosomes

sex_chromosome_DEGs <- function(all_genes,
                                ensemblID_colname="name"){

    # check input parameter is fine
    if(ensemblID_colname!="name"){
        if(!is.character(ensemblID_colname)){
            stop("Error: ensemblID_colname should be a string specifying the column name for the column containing ensembl IDs.")
        }
    }
    if(class(all_genes)!="list"){
        stop("Error: all_genes should be a list")
    }
    column_exists <- sapply(all_genes, function(x) ensemblID_colname %in% names(x))
    if(!all(column_exists)){
        stop("Error: every DGE analysis output in all_genes should contain the column with the name specified by ensemblID_colname")
    }

    # combine data for all celltypes into one dataframe
    combn_genes <- rbindlist(all_genes,idcol="celltype")
    # add symbol
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    genes <- unique(as.character(combn_genes$name))
    gene_IDs <- getBM(filters = "ensembl_gene_id", 
                      attributes = c("ensembl_gene_id","hgnc_symbol","chromosome_name"),
                      values = genes, mart = mart, useCache = FALSE)
    gene_IDs <- as.data.table(gene_IDs)
    setnames(gene_IDs,"ensembl_gene_id",ensemblID_colname)
    
    # remove any duplicates from reference set - two names for one ENSEMBL ID
    gene_IDs <- unique(gene_IDs,by=ensemblID_colname)
    setkey(combn_genes,name)
    # append gene names
    combn_genes[, gene_name := gene_IDs[combn_genes, on=.(name), x.hgnc_symbol]]
    # append chromosome
    combn_genes[, chromosome := gene_IDs[combn_genes, on=.(name), x.chromosome_name]]
    setorder(combn_genes,celltype,adj_pval,logFC)
    # include only genes on the sex chromosomes
    sexgenes <- combn_genes[combn_genes$chromosome == "X" | combn_genes$chromosome == "Y"]
    # remove unnecessary columns and reshape this as appropriate
    sexgenes <- split(sexgenes,sexgenes$celltype)
    sexgenes <- sexgenes[[names(sexgenes)]][, c("logFC","logCPM","LR","PValue","adj_pvalue","name")]

    # output DGE analysis subset to sex chromosome genes
    return(sexgenes)

}