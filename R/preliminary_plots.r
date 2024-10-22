# Define global variables
utils::globalVariables(c("..density.."))

#' Create preliminary plots for data exploration

#' @importFrom stats as.formula
#' @importFrom ggplot2 ggplot labs theme ggsave aes element_text element_blank geom_histogram
#' @importFrom cowplot theme_cowplot
#' @importFrom gridExtra arrangeGrob
#' @importFrom SingleCellExperiment colData

#' @param data the input data (should be an SCE object)
#' @param output_path base path in which outputs will be stored
#' @param sampleID sample ID
#' @param design the design formula of class type `formula`. Equation used to fit the model- data for the generalised linear model e.g. expression ~ sex + pmi + disease
#' @param sexID sex ID
#' @param celltypeID cell type ID
#' @param coeff which coefficient to carry out DE analysis with respect to
#' @param fdr the cut-off False Discovery Rate below which to select DEGs
#' @param y the column name in the SCE object for the return variable e.g. "diagnosis" - Case or disease. Default is the last variable in the design formula. y can be discrete (logistic regression) or continuous (linear regression)
#' @param region the column name in the SCE object for the study region. If there are multiple regions in the study (for example two brain regions). Pseudobulk values can be derived separately. Default is "single_region" which will not split by region.
#' @param control character specifying which control level for the differential expression analysis e.g. in a case/control/other study use "control" in the y column to compare against. NOTE only need to specify if more than two groups in y, leave as default value for two groups or continuous y. Default is NULL.
#' @param pval_adjust_method the adjustment method for the p-value in the differential expression analysis. Default is benjamini hochberg "BH". See  stats::p.adjust for available options
#' @param rmv_zero_count_genes whether genes with no count values in any cell should be removed. Default is TRUE

#' Saves all plots in the appropriate directory

preliminary_plots <- function(data,
                              output_path=getwd(),
                              sampleID="donor_id",
                              design="placeholder",
                              sexID="sex",
                              celltypeID="cell_type",
                              coeff="male",
                              fdr=0.05,
                              y=NULL,
                              region="single_region",
                              control=NULL,
                              pval_adjust_method="BH",
                              rmv_zero_count_genes=TRUE){

    # create output path if doesn't already exist
    setwd(output_path)
    dir.create(output_path,showWarnings=FALSE)
    # alter design
    if(design=="placeholder"){
        design=as.formula(paste0("~",sexID))
    }
    # validate function input params
    validate_input_parameters_power(data=data, output_path=output_path, sampleID=sampleID,
                                    design=design, sexID=sexID, celltypeID=celltypeID, 
                                    coeff=coeff, fdr=fdr, y=y,
                                    region=region, control=control, pval_adjust_method=pval_adjust_method,
                                    rmv_zero_count_genes=rmv_zero_count_genes)

    # get celltype name from dataset
    celltype_name <- toString(unique(data[[celltypeID]]))
    # check if DE analysis output present already in output_path
    if(!"DEout.RData" %in% list.files(output_path)){
    # run and save DE analysis
        assign("DEout", DGE_analysis(data, design=design, pseudobulk_ID=sampleID, celltype_ID=celltypeID, y=y, region=region, control=control, pval_adjust_method=pval_adjust_method, rmv_zero_count_genes=rmv_zero_count_genes, verbose=T, coef=coeff))
        save(DEout,file=paste0(output_path,"/DEout.RData"))
    }else{
        load(paste0(output_path,"/DEout.RData"))
    }
    
    ## plot histogram of effect sizes in DEGs and all genes
    # for DEGs (at fdr cut-off)
    DEGs_full <- subset(DEout$celltype_all_genes[[celltype_name]], adj_pval<fdr)
    logFCdist.plot <- ggplot(DEGs_full,aes(x=logFC)) + 
                            geom_histogram(aes(y=..density..), color="grey1", fill="grey1", alpha=0.5, bins=25) + 
                            labs(title=paste0("Distribution of LogFC (DEGs at ", fdr,"%)"),y="Density",x ="LogFC") + 
                            theme_cowplot() +
                            theme(axis.text=element_text(size=20),axis.title=element_text(size=22),plot.title=element_text(size=22))
    ggsave(paste0(output_path, "/logFCdist.png"),logFCdist.plot,width=20,height=20,units="cm",bg="white")
    ggsave(paste0(output_path, "/logFCdist.pdf"),logFCdist.plot,width=20,height=20,units="cm",bg="white")
    # for all genes
    allgenes <- DEout$celltype_all_genes[[celltype_name]]
    logFCdist_allgenes.plot <- ggplot(allgenes,aes(x=logFC)) + 
                            geom_histogram(aes(y=..density..), color="grey1", fill="grey1", alpha=0.5, bins=25) + 
                            labs(title=paste0("Distribution of LogFC (all genes)"),y="Density",x ="LogFC") + 
                            theme_cowplot() +
                            theme(axis.text=element_text(size=20),axis.title=element_text(size=22),plot.title=element_text(size=22))
    ggsave(paste0(output_path, "/logFCdist_allgenes.png"),logFCdist_allgenes.plot,width=20,height=20,units="cm",bg="white")
    ggsave(paste0(output_path, "/logFCdist_allgenes.pdf"),logFCdist_allgenes.plot,width=20,height=20,units="cm",bg="white")

    ## plot distribution of cells across individuals
    ## get table of sample ID/number of cells
    # get list of sample IDs
    coldata <- colData(data)
    IDs <- unique(coldata[[sampleID]])
    # get list of number of cells
    numIDs <- length(IDs)
    numCells <- list()
    for(i in 1:numIDs){
        numCells[[i]] <- dim(coldata[coldata[[sampleID]]==IDs[[i]],])[[1]]
    }
    numCells <- unlist(numCells)
    # put in dataframe
    cellcount <- data.frame(IDs)
    cellcount$numCells <- numCells
    # create and save histogram
    cellcount.plot <- ggplot(cellcount,aes(x=numCells)) + 
                            geom_histogram(aes(y=..density..), color="grey1", fill="grey1", alpha=0.5, bins=25) + 
                            labs(title="Distribution of numbers of cells across individuals",y="Density",x ="Number of cells") + 
                            theme_cowplot()+
                            theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=22),plot.title=element_text(size=22),axis.text.x = element_text(hjust=0.8))
    ggsave(paste0(output_path, "/cellcount.png"),cellcount.plot,width=20,height=20,units="cm",bg="white")
    ggsave(paste0(output_path, "/cellcount.pdf"),cellcount.plot,width=20,height=20,units="cm",bg="white")

    # put plots in a subplot together
    ggsave(paste0(output_path, "/QC_plots.png"),arrangeGrob(logFCdist.plot,logFCdist_allgenes.plot,cellcount.plot,ncol=3),width=60,height=20,units="cm",bg="white")
    ggsave(paste0(output_path, "/QC_plots.pdf"),arrangeGrob(logFCdist.plot,logFCdist_allgenes.plot,cellcount.plot,ncol=3),width=60,height=20,units="cm",bg="white")

}