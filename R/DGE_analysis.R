#'  Perform differential expression analysis using a pseudobulk approach on single-cell data
#'
#'  Conducts differential gene expression (DGE) analysis across cell types by aggregating single-cell RNA-seq data at the sample level (pseudobulk). This approach accounts for biological replicates and covariates, enabling robust identification of cell type–specific DEGs using generalized linear models. Returns a list of differential expression results for each cell type.

#' @importFrom data.table rbindlist
#' @importFrom qs qread
#' @importFrom stats as.formula

#' @param SCE A `SingleCellExperiment` object containing the input scRNA-seq data. You may also provide a path to an `.R`, `.rds`, or `.qs` file. If using a file, ensure the `SCE` object inside is named `SCE`.
#' @param design  A model formula specifying covariates for differential expression analysis. It should be of class `formula` (e.g., `~ sex + pmi + disease`). This formula is used to fit a generalized linear model.
#' @param sampleID Name of the column in the `SCE` metadata that identifies biological replicates (e.g., patient ID). This column is used for grouping in the pseudobulk approach.
#' @param celltypeID Name of the column in the `SCE` metadata indicating cell type labels. This is used to identify celltype specific DEGs. If there is only one cell type in the analysis, the analysis will automatically adjust.
#' @param y Name of the column in the `SCE` metadata representing the response variable (e.g., "diagnosis" - case or disease). If not specified, defaults to the last variable in the `design` formula. Accepts both categorical (logistic regression) and continuous (linear regression) variables.
#' @param coef Character string indicating the level of the response variable (`y`) to test for in differential expression. For case-control studies, this would typically be "case" (e.g. "AD"). Typically used in binary comparisons. Not required for continuous outcomes.
#' @param region Optional column in `SCE` metadata indicating the tissue or brain region. If present, differential expression is performed within each region separately. Defaults to "single_region" (i.e., no regional split).
#' @param control  Optional. Character string specifying the control level in the response variable (`y`) to compare against. Only required if `y` contains more than two levels. Ignored for binary or continuous outcomes.
#' @param pval_adjust_method Method used to adjust p-values for multiple testing. Default is "BH" (Benjamini–Hochberg). See `stats::p.adjust` for available options.
#' @param adj_pval Adjusted p-value threshold (0–1) used to define significance in the differential expression results. Default is 0.05.
#' @param output_path  Directory where output plots from the DE analysis will be saved. If set to `FALSE`, no plots will be generated. Defaults to a folder named "sc.cell.type.de.graphs" in the working directory.
#' @param verbose Logical. Whether to print progress and additional messages during the differential analysis analysis. Default is `FALSE`.
#' @param rmv_zero_count_genes Logical. Whether to remove genes with zero counts across all cells. Default is `TRUE`.
#' @param save Whether to save DE analysis results as an `.RData` file. Default is `FALSE`. Note: power analysis and downsampling functions handle saving separately.

#' @return A list containing:
#'             celltype_DEGs: list with the differentially expressed genes (DEGs) for each cell type
#'             celltype_all_genes: list with all genes along with their differential expression scores for each cell type
#'             celltype_counts: vector with the counts of cells after QC in each cell type
#' @export

#' @examples
#'\dontrun{
#'# To use the function input the formula of the comparison you want to make along with the names
#'# of the pseudobulk ID and celltype ID
#'# If you want to compare disease and control across cell types, taking into account sex use:
#'#  formula = ~ sex + disease. Firstly load  your SCE object (can otherwise specify location)
#' library(qs)
#' library(SingleCellExperiment)
#' SCE <- qs::qread("sce.qs")
#' DGE_analysis.return_incl_sex<- DGE_analysis(SCE_small,
#' design= ~ sex + pathological_diagnosis_original,
#' sampleID="sample_id", celltypeID="allan_celltype",coef="AD")
#'# If you only want to also account for other variables as well as sex, such as
#'# postmortem interval (PMI):
#' DGE_analysis.return_sex_pmi<- DGE_analysis(SCE_small,
#' design= ~ sex + PMI + pathological_diagnosis_original,
#' sampleID="sample_id", celltypeID="allan_celltype",coef="AD")
#' # If you only want to compare disease and control across cell types w/o any other constraints:
#' DGE_analysis.return<- DGE_analysis(SCE_small,design= ~ pathological_diagnosis_original,
#' sampleID="sample_id", celltypeID="allan_celltype",coef="AD")
#' # A good way to validate the DEG analysis approach is to run a sex comparison
#' # You would expect genes on the sex chromosomes to be DEGs, this can be done:
#' DGE_analysis.sex.return<- DGE_analysis(SCE,design= ~ sex,
#' sampleID="sample_id", celltypeID="allan_celltype",coef="M")
#' #get gene names of chromosomal DEGs
#' library('biomaRt')
#' mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#' genes <- unique(as.character(our_degs$name))
#' gene_IDs <- getBM(filters= "ensembl_gene_id",
#'                   attributes= c("ensembl_gene_id","chromosome_name","hgnc_symbol"),
#'                   values = genes, mart= mart,useCache = FALSE)
#' gene_IDs <- as.data.table(gene_IDs)
#' setnames(gene_IDs,"ensembl_gene_id","name")
#' DEGs <- data.table::rbindlist(sc.cell.type.de.return$celltype_DEGs,idcol="celltype")
#' setkey(DEGs,name)
#' #append gene names
#' DEGs[, gene_name := gene_IDs[DEGs, on=.(name), x.hgnc_symbol]]
#' DEGs[, chromosome_name := gene_IDs[DEGs, on=.(name), x.chromosome_name]]
#' #Xist should be there, key indicator of correct DEG analysis
#' DEGs[gene_name=="XIST",c("celltype","adj_pval","lfc")]
#' #Present in all cell types
#' # There should also be high number of other DEGs on sex chromosomes, check:
#' our_degs[chromosome_name %in% c("X","Y"),
#'          .N,by=.(celltype,chromosome_name)]
#'}
#'
#' # Simple runnable example with bundled data
#' micro_tsai <- system.file("extdata", "Tsai_Micro.qs", package="poweranalysis")
#' SCE_tsai <- qs::qread(micro_tsai)
#' DGE_tsai <- poweranalysis::DGE_analysis(
#'     SCE_tsai,
#'     design = ~ sex,
#'     celltypeID = "cluster_celltype",
#'     sampleID = "sample_id",
#'     output_path = tempdir(),
#'     coef = "M"
#' )
#' DGE_tsai

DGE_analysis <- function(SCE, design, sampleID, celltypeID, y=NULL,
                         region="single_region", coef=NULL, control=NULL,
                         pval_adjust_method = "BH", adj_pval=0.05,
                         output_path="sc.cell.type.de.graphs/",
                         rmv_zero_count_genes=T, verbose=F, save=F){

    # need to load SCE if a directory is passed
    if(class(SCE)[1]=="character"){
        if(!file.exists(SCE))
            stop("Directory to SCE file doesn't exist, please check your file location")
        # load the dataset
        if(substr(SCE,nchar(SCE)-2,nchar(SCE))==".qs"){
            # load QS objects
            SCE <- qread(SCE)
        }else if (substr(SCE,nchar(SCE)-3,nchar(SCE))==".rds"){
            # rds object
            SCE <- readRDS(SCE)
        }else{
            #normal R object - saved object must also be called SCE for this to work, will hit error in input check otherwise
            load(SCE)
        }
    }

    # sense check inputs
    validate_input_parameters_de(SCE, design, sampleID, celltypeID, y,
                                 region, coef, control, pval_adjust_method,
                                 adj_pval, output_path, rmv_zero_count_genes,
                                 verbose)

    # get counts of each cell type
    #counts(SCE) <- as.matrix(counts(SCE))
    counts_celltypes <- SCE[[celltypeID]]
    counts_celltypes <- as.vector(table(counts_celltypes))
    names(counts_celltypes) <- names(table(SCE[[celltypeID]]))

    # first format formula
    design_txt <- paste0(deparse(design,width.cutoff = 500),collapse=',')
    # make design formula minus anything before ~
    formula <- as.formula(gsub(".*~","~",design_txt))
    # if exists remove everything before ~ this format is necessary for glm_gp()
    design_txt <- gsub(".*~","",design_txt)
    # split design by "+" to get components to bring forward for pseudobulk, add in pseudobulk by column
    pb_columns <- c(gsub("^\\s+|\\s+$","",strsplit(design_txt, "[+]")[[1]]),
                        sampleID)
    # if y not specified take last value in design matrix
    if(is.null(y))
        y <- design_txt[[length(design_txt)]]
    # check if continuous or categorical variable to be modeled
    y_contin <- FALSE
    if(is.numeric(SCE[[y]]))
        y_contin <- TRUE

    # get pseudobulk values
    celltypes <- unique(SCE[[celltypeID]])
    if(isTRUE(verbose))
        message("Deriving pseudobulk data")
    pb_dat <-
        lapply(celltypes,function(x)
            make_pseudobulk(SCE[,SCE[[celltypeID]]==x],
                            sampleID=sampleID,
                            pb_columns=pb_columns,
                            region=region,
                            rmv_zero_count_genes=rmv_zero_count_genes))
    names(pb_dat) <- celltypes

    # run edgeR LRT DE analysis
    celltype_de <-
        differential_expression(pb_dat,formula,y_name=y,y_contin,coef,control,
                                pval_adjust_method, adj_pval, verbose)
    # get sig DEGs for each
    celltype_DEGs <- lapply(celltype_de, function(x) x[x$adj_pval<adj_pval,])

    unique_genes <-lapply(celltype_DEGs, function(x) x$name)
    unique_degs <- unique(unlist(unique_genes))

    # if no DEGs found break but return the DE analysis scores still
    if(length(unique_degs)==0)
        return(list("celltype_DEGs"=celltype_DEGs,
                    "celltype_all_genes"=celltype_de,
                    "celltype_counts"=counts_celltypes))

    if(isTRUE(verbose))
        message(length(unique_degs)," unique DEGs foundacross all cell types")

    celltype_all_genes_dt <- data.table::rbindlist(celltype_de,idcol = T)
    setnames(celltype_all_genes_dt,".id","celltype")
    celltype_DEGs_dt <- data.table::rbindlist(celltype_DEGs,idcol = T)
    setnames(celltype_DEGs_dt,".id","celltype")

    if(!isFALSE(output_path)){
        if(isTRUE(verbose))
        message("Plotting the results of differential expression analysis")
        #make plots for DE analysis
        plot_de_analysis(pb_dat,y,celltype_DEGs_dt,celltype_all_genes_dt,
                         counts_celltypes,output_path)
        # save the DE analysis results if specified
        DEout <- list("celltype_DEGs"=celltype_DEGs,
                      "celltype_all_genes"=celltype_de,
                      "celltype_counts"=counts_celltypes)
        if(isTRUE(save))
            save(DEout, file = file.path(output_path,"DEout.RData"))
    }

    # return the DE analysis results
    return(DEout)

}
