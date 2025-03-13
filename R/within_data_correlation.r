# Define global variables
utils::globalVariables(c("DEout"))

#' Create correlation plots of the effect sizes, between random permutations and subsets of a given dataset

#' @importFrom ggplot2 ggsave

#' @param SCE the input data (should be an SCE object)
#' @param output_path base path in which outputs will be stored
#' @param sampleID sample ID
#' @param design the design formula of class type `formula`. Equation used to fit the model- data for the generalised linear model e.g. expression ~ sex + pmi + disease
#' @param sexID sex ID
#' @param celltypeID cell type ID
#' @param coeff which coefficient to carry out DE analysis with respect to
#' @param N_randperms number of randomised permutations of the dataset (based on sex) to be correlated
#' @param N_subsetpairs number of pairs of subsets of the dataset to be correlated
#' @param y the column name in the SCE object for the return variable e.g. "diagnosis" - Case or disease. Default is the last variable in the design formula. y can be discrete (logistic regression) or continuous (linear regression)
#' @param region the column name in the SCE object for the study region. If there are multiple regions in the study (for example two brain regions). Pseudobulk values can be derived separately. Default is "single_region" which will not split by region.
#' @param control character specifying which control level for the differential expression analysis e.g. in a case/control/other study use "control" in the y column to compare against. NOTE only need to specify if more than two groups in y, leave as default value for two groups or continuous y. Default is NULL.
#' @param pval_adjust_method the adjustment method for the p-value in the differential expression analysis. Default is benjamini hochberg "BH". See stats::p.adjust for available options
#' @param rmv_zero_count_genes whether genes with no count values in any cell should be removed. Default is TRUE

#' Saves all plots in the appropriate directory

#' @export

within_data_correlation <- function(SCE,
                                    output_path=getwd(),
                                    sampleID="donor_id",
                                    design="placeholder",
                                    sexID="sex",
                                    celltypeID="cell_type",
                                    coeff="male",
                                    N_randperms=5,
                                    N_subsetpairs=5,
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
    validate_input_parameters_power(SCE=SCE, output_path=output_path, sampleID=sampleID,
                                    design=design, sexID=sexID, celltypeID=celltypeID,
                                    coeff=coeff, N_randperms=N_randperms, N_subsetpairs=N_subsetpairs,
                                    y=y, region=region, control=control,
                                    pval_adjust_method=pval_adjust_method, rmv_zero_count_genes=rmv_zero_count_genes)

    # check if DE analysis output present already in output_path
    if(!"DEout.RData" %in% list.files(output_path)){
        stop("Error: DGE analysis output file (DEout.RData) for the full dataset is not in the specified (output_path) directory!")
    }else{
        load(file.path(output_path,"DEout.RData"))
    }
    # create directory for correlation analysis if doesn't already exist
    dir.create(file.path(output_path, "corr_analysis"), showWarnings=FALSE)

    # get all genes from main dataset
    allgenes_full <- DEout$celltype_all_genes
    # list of celltypes, studies
    celltypes <- names(allgenes_full)
    ## create random permutations and subsets of data, run DE analysis on them
    # random perms
    permutedData <- random_permutations(SCE,sampleID,sexID,N_randperms)
    # list of data
    allstudies <- list()
    names <- c()
    # run DGE analysis on these
    for(perm in 1:N_randperms){
        # create relevant directory and move to it
        dir.create(file.path(output_path, paste0("corr_analysis/randperm",toString(perm))),showWarnings=FALSE)
        savepath <- file.path(output_path, paste0("corr_analysis/randperm",toString(perm)))
        # run DE analysis
        assign(paste0("perm_",toString(perm),"DEout"), DGE_analysis(permutedData[[perm]], design=design, sampleID=sampleID, celltypeID=celltypeID, y=y, region=region, control=control, pval_adjust_method=pval_adjust_method, rmv_zero_count_genes=rmv_zero_count_genes, verbose=T, coef=coeff))
        # save
        save(list=eval(paste0("perm_",toString(perm),"DEout")),file=file.path(savepath,paste0("perm_",toString(perm),"DEout.RData")))
        # get all genes
        assign(paste0("rand",toString(perm),"_genes"),eval(as.name(paste0("perm_",toString(perm),"DEout")))$celltype_all_genes)
        # add to list
        allstudies <- append(allstudies, list(eval(as.name(paste0("rand",toString(perm),"_genes")))))
        # add name
        names <- c(names, paste0("rand_perm_",toString(perm)))
    }
    # add main study
    allstudies <- append(allstudies, list(allgenes_full))
    # add name
    names <- c(names, "full_dataset")
    # subsets
    subsets <- subset_pairs(SCE,sampleID,N_subsetpairs)
    for(subset in 1:N_subsetpairs){
        ## subset a
        # create relevant directory and move to it
        dir.create(file.path(output_path, paste0("corr_analysis/subset",toString(subset),"a")),showWarnings=FALSE)
        savepath_a <- file.path(output_path, paste0("corr_analysis/subset",toString(subset),"a"))
        # run DE analysis
        assign(paste0("subset",toString(subset),"a_DEout"), DGE_analysis(subsets[[subset]][[1]], design=design, sampleID=sampleID, celltypeID=celltypeID, y=y, region=region, control=control, pval_adjust_method=pval_adjust_method, rmv_zero_count_genes=rmv_zero_count_genes, verbose=T, coef=coeff))
        # save
        save(list=eval(paste0("subset",toString(subset),"a_DEout")),file=file.path(savepath_a,paste0("subset",toString(subset),"a_DEout.RData")))
        ## subset b
        # create relevant directory and move to it
        dir.create(file.path(output_path, paste0("corr_analysis/subset",toString(subset),"b")),showWarnings=FALSE)
        savepath_b <- file.path(output_path, paste0("corr_analysis/subset",toString(subset),"b"))
        # run DE analysis
        assign(paste0("subset",toString(subset),"b_DEout"), DGE_analysis(subsets[[subset]][[2]], design=design, sampleID=sampleID, celltypeID=celltypeID, y=y, region=region, control=control, pval_adjust_method=pval_adjust_method, rmv_zero_count_genes=rmv_zero_count_genes, verbose=T, coef=coeff))
        # save
        save(list=eval(paste0("subset",toString(subset),"b_DEout")),file=file.path(savepath_b,paste0("subset",toString(subset),"b_DEout.RData")))
        # get all genes
        assign(paste0("subset",toString(subset),"a_genes"),eval(as.name(paste0("subset",toString(subset),"a_DEout")))$celltype_all_genes)
        assign(paste0("subset",toString(subset),"b_genes"),eval(as.name(paste0("subset",toString(subset),"b_DEout")))$celltype_all_genes)
        allstudies <- append(allstudies, list(eval(as.name(paste0("subset",toString(subset),"a_genes")))))
        allstudies <- append(allstudies, list(eval(as.name(paste0("subset",toString(subset),"b_genes")))))
        # add name
        names <- c(names, paste0("subset",toString(subset),"a"))
        names <- c(names, paste0("subset",toString(subset),"b"))
    }

    ## correlation analysis
    # remove non-overlapping celltypes
    names(allstudies) <- names
    allstudies <- lapply(allstudies,"[",celltypes)
    ## correlations
    # 0.05
    savepath_corrs <- file.path(output_path, "corr_analysis")
    corr_0.05pval <- plot_mean_correlation("full_dataset",allstudies,celltypes,0.05)
    ggsave(file.path(savepath_corrs,"corrPlot_0.05pval.png"),corr_0.05pval[[1]],width=20,height=20,units="cm",bg="white")
    ggsave(file.path(savepath_corrs,"corrPlot_0.05pval.pdf"),corr_0.05pval[[1]],width=20,height=20,units="cm",bg="white")
    # 0.025
    corr_0.025pval <- plot_mean_correlation("full_dataset",allstudies,celltypes,0.025)
    ggsave(file.path(savepath_corrs,"corrPlot_0.025pval.png"),corr_0.025pval[[1]],width=20,height=20,units="cm",bg="white")
    ggsave(file.path(savepath_corrs,"corrPlot_0.025pval.pdf"),corr_0.025pval[[1]],width=20,height=20,units="cm",bg="white")
    # 0.01
    corr_0.01pval <- plot_mean_correlation("full_dataset",allstudies,celltypes,0.01)
    ggsave(file.path(savepath_corrs,"corrPlot_0.01pval.png"),corr_0.01pval[[1]],width=20,height=20,units="cm",bg="white")
    ggsave(file.path(savepath_corrs,"corrPlot_0.01pval.pdf"),corr_0.01pval[[1]],width=20,height=20,units="cm",bg="white")
    # 0.001
    corr_0.001pval <- plot_mean_correlation("full_dataset",allstudies,celltypes,0.001)
    ggsave(file.path(savepath_corrs,"corrPlot_0.001pval.png"),corr_0.001pval[[1]],width=20,height=20,units="cm",bg="white")
    ggsave(file.path(savepath_corrs,"corrPlot_0.001pval.pdf"),corr_0.001pval[[1]],width=20,height=20,units="cm",bg="white")
    # 0.0001
    corr_0.0001pval <- plot_mean_correlation("full_dataset",allstudies,celltypes,0.0001)
    ggsave(file.path(savepath_corrs,"corrPlot_0.0001pval.png"),corr_0.0001pval[[1]],width=20,height=20,units="cm",bg="white")
    ggsave(file.path(savepath_corrs,"corrPlot_0.0001pval.pdf"),corr_0.0001pval[[1]],width=20,height=20,units="cm",bg="white")
    
}