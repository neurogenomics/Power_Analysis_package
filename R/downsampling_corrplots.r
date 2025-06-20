#' Create correlation plots of the effect sizes of the top 1000 and top 500 DEGs

#' @importFrom gtools mixedsort
#' @importFrom stringr str_sub
#' @importFrom ggplot2 labs theme ggsave element_text scale_fill_gradient2
#' @importFrom ggcorrplot ggcorrplot

#' @param SCE the input data (should be an SCE object)
#' @param range_downsampled vector or list containing values which the data will be downsampled at, in ascending order
#' @param output_path base path in which outputs will be stored
#' @param inpath base path where downsampled DGE analysis output folders are stored (taken to be output_path if not provided)
#' @param sampled downsampling carried out based on what (either "individuals" or "cells")
#' @param sampleID sample ID
#' @param celltypeID cell type ID
#' @param coeff which coefficient to carry out DE analysis with respect to
#' @param Nperms number of subsets created when downsampling at each level
#' @param y the column name in the SCE object for the return variable e.g. "diagnosis" - Case or disease. Default is the last variable in the design formula. y can be discrete (logistic regression) or continuous (linear regression)
#' @param region the column name in the SCE object for the study region. If there are multiple regions in the study (for example two brain regions). Pseudobulk values can be derived separately. Default is "single_region" which will not split by region.
#' @param control character specifying which control level for the differential expression analysis e.g. in a case/control/other study use "control" in the y column to compare against. NOTE only need to specify if more than two groups in y, leave as default value for two groups or continuous y. Default is NULL.
#' @param pval_adjust_method the adjustment method for the p-value in the differential expression analysis. Default is benjamini hochberg "BH". See  stats::p.adjust for available options
#' @param rmv_zero_count_genes whether genes with no count values in any cell should be removed. Default is TRUE

#' Saves all plots in the appropriate directory

downsampling_corrplots <- function(SCE,
                                   range_downsampled="placeholder",
                                   output_path=getwd(),
                                   inpath="placeholder",
                                   sampled="individuals",
                                   sampleID="donor_id",
                                   celltypeID="cell_type",
                                   coeff="male",
                                   Nperms=20,
                                   y=NULL,
                                   region="single_region",
                                   control=NULL,
                                   pval_adjust_method="BH",
                                   rmv_zero_count_genes=TRUE){
    
    # alter range_downsampled
    if(identical(range_downsampled,"placeholder")){
        range_downsampled <- downsampling_range(SCE, sampled, sampleID)
    }
    # alter inpath
    if(inpath=="placeholder"){
        inpath <- output_path
    }

    # create output path if doesn't already exist
    setwd(output_path)
    dir.create(output_path,showWarnings=FALSE,recursive=TRUE)
    # create directory for correlation analysis
    dir.create(file.path(output_path, "corr_analysis"), showWarnings=FALSE,recursive=TRUE)
    
    # get celltype name from dataset
    celltype_name <- toString(unique(SCE[[celltypeID]]))
    # check if DE analysis output present already in output_path
    if(!"DEout.RData" %in% list.files(output_path)){
        # run and save DE analysis
        assign("DEout", DGE_analysis(SCE, design=design, sampleID=sampleID, celltypeID=celltypeID, y=y, region=region, control=control, pval_adjust_method=pval_adjust_method, rmv_zero_count_genes=rmv_zero_count_genes, verbose=T, coef=coeff))
        save(DEout,file=file.path(output_path,"DEout.RData"))
    }else{
        load(file.path(output_path,"DEout.RData"))
    }

    ## correlation using down-sampled datasets - top 1000 genes
    allgenes <- DEout$celltype_all_genes[[celltype_name]]
    # rank DEGs by significance
    allgenes_ordered <- allgenes[order(allgenes$PValue),]
    # get top 1000 most signficant DEGs from full dataset
    top1000_full <- allgenes_ordered[1:1000,]
    top1000 <- data.frame(top1000_full$name)
    colnames(top1000) <- c("name")
    top1000$`logFC` <- top1000_full$logFC
    ## correlation using down-sampled datasets - top 500 genes
    # get top 500 most signficant DEGs from full dataset
    top500_full <- allgenes_ordered[1:500,]
    top500 <- data.frame(top500_full$name)
    colnames(top500) <- c("name")
    top500$`logFC` <- top500_full$logFC

    ## if down-sampled individuals
    if(sampled=="individuals"){
        
        # define path and folders
        path <- file.path(inpath,"DE_downsampling")
        downsampled_folders <- paste0(paste(range_downsampled, sep=" ", collapse=NULL),"samples")
        # get corr matrices for each permutation
        corrMats_1000 <- list()
        for(i in 1:Nperms){
            setwd(path)
            for(folder in mixedsort(list.files())){
                # remove other files above
                if(folder%in%downsampled_folders){
                    newpath <- file.path(path,folder)
                    # enter directory
                    setwd(newpath)
                    # get num_samples
                    num_samples <- str_sub(folder,1,-8) 
                    # go into permutation
                    subpath <- file.path(newpath,paste0(num_samples,"_",i))
                    setwd(subpath)
                    # read DGE analysis output
                    load(paste0("DEout",num_samples,"_",i,".RData"))
                    # get df for top 1000 genes
                    all_genes <- get(paste0("DEout_",num_samples))$celltype_all_genes[[celltype_name]]
                    assign(paste0("samples_",num_samples),subset(all_genes,name%in%top1000$name))
                    # move out again
                    setwd(newpath)
                }
            }
            allstudies <- mget(paste0("samples_", paste(range_downsampled, sep=" ", collapse=NULL)))
            names(allstudies) <- paste0(paste(range_downsampled, sep=" ", collapse=NULL),"_samples")
            # compute correlation for this permutation
            corrMats_1000[[i]] <- compute_downsampled_corr(allstudies)[[1]]
        }
        # compute mean of correlation matrices
        meanCorr_1000 <- Reduce("+",corrMats_1000) / length(corrMats_1000)
        # plot
        setwd(file.path(output_path, "corr_analysis"))
        # plot correlation matrix
        meanCorr_downsampling_1000.plot <- ggcorrplot(round(meanCorr_1000,2), 
                hc.order = F, insig="pch",pch=5,pch.col = "grey",
                pch.cex=9,
                title=paste0("Mean LogFC Correlation Matrix (across permutations) for top 1000 genes when downsampling"),
                colors = c("#FC4E07", "white", "#00AFBB"),
                outline.color = "white", lab = TRUE, lab_size=3.5,
                sig.level=0.05) +
                theme(plot.title = element_text(hjust = 0.7)) +
                scale_fill_gradient2(breaks=c(0, 1), limit=c(0, 1)) +
                labs(fill="Corr")
        ggsave("meanCorr_downsampling_1000.png",meanCorr_downsampling_1000.plot,width=25,height=20,units="cm",bg="white")
        ggsave("meanCorr_downsampling_1000.pdf",meanCorr_downsampling_1000.plot,width=25,height=20,units="cm",bg="white")

        # get corr matrices for each permutation
        corrMats_500 <- list()
        for(i in 1:Nperms){
            setwd(path)
            for(folder in mixedsort(list.files())){
                # remove other files above
                if(folder%in%downsampled_folders){
                    newpath <- file.path(path,folder)
                    # enter directory
                    setwd(newpath)
                    # get num_samples
                    num_samples <- str_sub(folder,1,-8) 
                    # go into permutation
                    subpath <- file.path(newpath,paste0(num_samples,"_",i))
                    setwd(subpath)
                    # read DGE analysis output
                    load(paste0("DEout",num_samples,"_",i,".RData"))
                    # get df for top 500 genes
                    all_genes <- get(paste0("DEout_",num_samples))$celltype_all_genes[[celltype_name]]
                    assign(paste0("samples_",num_samples),subset(all_genes,name%in%top500$name))
                    # move out again
                    setwd(newpath)
                }
            }
            allstudies <- mget(paste0("samples_", paste(range_downsampled, sep=" ", collapse=NULL)))
            names(allstudies) <- paste0(paste(range_downsampled, sep=" ", collapse=NULL),"_samples")
            # compute correlation for this permutation
            corrMats_500[[i]] <- compute_downsampled_corr(allstudies)[[1]]
        }
        # compute mean of correlation matrices
        meanCorr_500 <- Reduce("+",corrMats_500) / length(corrMats_500)
        # plot
        setwd(file.path(output_path, "corr_analysis"))
        # plot correlation matrix
        meanCorr_downsampling_500.plot <- ggcorrplot(round(meanCorr_500,2), 
                hc.order = F, insig="pch",pch=5,pch.col = "grey",
                pch.cex=9,
                title=paste0("Mean LogFC Correlation Matrix (across permutations) for top 500 genes when downsampling"),
                colors = c("#FC4E07", "white", "#00AFBB"),
                outline.color = "white", lab = TRUE, lab_size=3.5,
                sig.level=0.05) +
                theme(plot.title = element_text(hjust = 0.7)) +
                scale_fill_gradient2(breaks=c(0, 1), limit=c(0, 1)) +
                labs(fill="Corr")
        ggsave("meanCorr_downsampling_500.png",meanCorr_downsampling_500.plot,width=25,height=20,units="cm",bg="white")
        ggsave("meanCorr_downsampling_500.pdf",meanCorr_downsampling_500.plot,width=25,height=20,units="cm",bg="white")

    }else{

        ## if down-sampled cells
        # define path and folders
        path <- file.path(inpath,"DE_downsampling_cells")
        downsampled_folders <- paste0(paste(range_downsampled, sep=" ", collapse=NULL),"cells_persample")
        # get corr matrices for each permutation
        corrMats_cells_1000 <- list()
        for(i in 1:Nperms){
            setwd(path)
            for(folder in mixedsort(list.files())){
                # remove other files above
                if(folder%in%downsampled_folders){
                    newpath <- file.path(path,folder)
                    # enter directory
                    setwd(newpath)
                    # get num_cells
                    num_cells <- str_sub(folder,1,-16)
                    # go into permutation
                    subpath <- file.path(newpath,paste0(num_cells,"_",i))
                    setwd(subpath)
                    # read DGE analysis output
                    load(paste0("DEout",num_cells,"_",i,".RData"))
                    # get df for top 1000 genes
                    all_genes <- get(paste0("DEout_",num_cells))$celltype_all_genes[[celltype_name]]
                    assign(paste0("cells_",num_cells),subset(all_genes,name%in%top1000$name))
                    # move out again
                    setwd(newpath)
                }
            }
            allstudies <- mget(paste0("cells_", paste(range_downsampled, sep=" ", collapse=NULL)))
            names(allstudies) <- paste0(paste(range_downsampled, sep=" ", collapse=NULL),"_cells/sample")
            # compute correlation for this permutation
            corrMats_cells_1000[[i]] <- compute_downsampled_corr(allstudies,"cells")[[1]]
        }
        # compute mean of correlation matrices
        meanCorr_cells_1000 <- Reduce("+",corrMats_cells_1000) / length(corrMats_cells_1000)
        # plot
        setwd(file.path(output_path, "corr_analysis"))
        # plot correlation matrix
        meanCorr_downsampling_cells_1000.plot <- ggcorrplot(round(meanCorr_cells_1000,2), 
                hc.order = F, insig="pch",pch=5,pch.col = "grey",
                pch.cex=9,
                title=paste0("Mean LogFC Correlation Matrix (across permutations) for top 1000 genes when downsampling cells"),
                colors = c("#FC4E07", "white", "#00AFBB"),
                outline.color = "white", lab = TRUE, lab_size=3.5,
                sig.level=0.05) +
                theme(plot.title = element_text(hjust = 0.7)) +
                scale_fill_gradient2(breaks=c(0, 1), limit=c(0, 1)) +
                labs(fill="Corr")
        ggsave("meanCorr_downsampling_cells_1000.png",meanCorr_downsampling_cells_1000.plot,width=25,height=20,units="cm",bg="white")
        ggsave("meanCorr_downsampling_cells_1000.pdf",meanCorr_downsampling_cells_1000.plot,width=25,height=20,units="cm",bg="white")

        ## correlation using down-sampled datasets (cells) - top 500 genes
        # load in DGE analysis outputs for down-sampled datasets, take union of DEGs across perms
        setwd(path)
        # get corr matrices for each permutation
        corrMats_cells_500 <- list()
        for(i in 1:Nperms){
            setwd(path)
            for(folder in mixedsort(list.files())){
                # remove other files above
                if(folder%in%downsampled_folders){
                    newpath <- file.path(path,folder)
                    # enter directory
                    setwd(newpath)
                    # get num_cells
                    num_cells <- str_sub(folder,1,-16)
                    # go into permutation
                    subpath <- file.path(newpath,paste0(num_cells,"_",i))
                    setwd(subpath)
                    # read DGE analysis output
                    load(paste0("DEout",num_cells,"_",i,".RData"))
                    # get df for top 500 genes
                    all_genes <- get(paste0("DEout_",num_cells))$celltype_all_genes[[celltype_name]]
                    assign(paste0("cells_",num_cells),subset(all_genes,name%in%top500$name))
                    # move out again
                    setwd(newpath)
                }
            }
            allstudies <- mget(paste0("cells_", paste(range_downsampled, sep=" ", collapse=NULL)))
            names(allstudies) <- paste0(paste(range_downsampled, sep=" ", collapse=NULL),"_cells/sample")
            # compute correlation for this permutation
            corrMats_cells_500[[i]] <- compute_downsampled_corr(allstudies,"cells")[[1]]
        }
        # compute mean of correlation matrices
        meanCorr_cells_500 <- Reduce("+",corrMats_cells_500) / length(corrMats_cells_500)
        # plot
        setwd(file.path(output_path, "corr_analysis"))
        # plot correlation matrix
        meanCorr_downsampling_cells_500.plot <- ggcorrplot(round(meanCorr_cells_500,2), 
                hc.order = F, insig="pch",pch=5,pch.col = "grey",
                pch.cex=9,
                title=paste0("Mean LogFC Correlation Matrix (across permutations) for top 500 genes when downsampling cells"),
                colors = c("#FC4E07", "white", "#00AFBB"),
                outline.color = "white", lab = TRUE, lab_size=3.5,
                sig.level=0.05) +
                theme(plot.title = element_text(hjust = 0.7)) +
                scale_fill_gradient2(breaks=c(0, 1), limit=c(0, 1)) +
                labs(fill="Corr")
        ggsave("meanCorr_downsampling_cells_500.png",meanCorr_downsampling_cells_500.plot,width=25,height=20,units="cm",bg="white")
        ggsave("meanCorr_downsampling_cells_500.pdf",meanCorr_downsampling_cells_500.plot,width=25,height=20,units="cm",bg="white")

    }
    
}