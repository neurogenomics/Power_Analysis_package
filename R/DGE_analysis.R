#' Perform differential expression analysis across cell types based on single cell data with Pseudobulk approach. Returns the results of this differential expression analysis

#' @importFrom data.table rbindlist
#' @importFrom qs qread
#' @importFrom stats as.formula

#' @param SCE SingleCellExperiment object, a specialised S4 class for storing data from single-cell experiments. A location of an R, rds or qs file to be loaded can also be passed. If using R objects make sure SCE object saved with the name SCE
#' @param design Design formula of class type `formula`. Equation used to fit the model- data for the generalised linear model e.g. expression ~ sex + pmi + disease.
#' @param sampleID Column name in the SCE object to perform pseudobulk on, usually the patient identifier. This column is used for grouping in the pseudobulk approach
#' @param celltypeID Column name in the SCE object for the cell type variable. This is used to identify celltype specific DEGs. If there is only one cell type in the analysis, the analysis will be automatically updated accordingly.
#' @param y Column name in the SCE object for the return variable e.g. "diagnosis" - Case or disease. Default is the last variable in the design formula. y can be discrete (logisitc regression) or continuous (linear regression)
#' @param region Column name in the SCE object for the study region. If there are multiple regions in the study (for example two brain regions). Pseudobulk values can be derived separately. Default is "single_region" which will not split by region.
#' @param coef Character specifying which level to investigate for the differential expression analysis e.g. in a case/control study use "case" if case is the identifier in the y column to get positive fold changes to relate to case samples. leave as default value for continuous y. Default is NULL.
#' @param control Character specifying which control level for the differential expression analysis e.g. in a case/control/other study use "control" in the y column to compare against. NOTE only need to specify if more than two groups in y, leave as default value for two groups or continuous y. Default is NULL.
#' @param pval_adjust_method Adjustment method for the p-value in the differential expression analysis. Default is benjamini hochberg "BH". See  stats::p.adjust for available options
#' @param adj_pval Adjusted p-value cut-off for the differential expression analysis, 0-1 range
#' @param output_path Folder where the graphs from the differential expression analysis are saved. Default will create a folder in the current working directory "sc.cell.type.de.graphs". False will skip plotting.
#' @param verbose Logical indicating if extra information about the differential expression analysis should be printed
#' @param rmv_zero_count_genes Whether genes with no count values in any cell should be removed. Default is TRUE
#' @param save Whether or not to save the DE analysis results as an RData file. Default is FALSE (down-sampling/power analysis saves these outputs itself)

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
