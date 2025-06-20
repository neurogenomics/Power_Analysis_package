# Define global variables
utils::globalVariables(c("design","PValue","logFC","name","variable"))

#' Create plots for power analysis, with down-sampling based either on the individuals or cells

#' @importFrom reshape2 melt
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot labs theme ggsave ggtitle aes element_text geom_boxplot
#' @importFrom cowplot theme_cowplot
#' @importFrom viridis scale_fill_viridis
#' @importFrom gtools mixedsort
#' @importFrom stringr str_sub
#' @importFrom gridExtra arrangeGrob

#' @param SCE the input data (should be an SCE object)
#' @param range_downsampled vector or list containing values which the data will be downsampled at, in ascending order
#' @param output_path base path in which outputs will be stored
#' @param inpath base path where downsampled DGE analysis output is stored (taken to be output_path if not provided)
#' @param sampled downsampling carried out based on what (either "individuals" or "cells")
#' @param sampleID sample ID
#' @param celltypeID cell type ID
#' @param coeff which coefficient to carry out DE analysis with respect to
#' @param fdr the cut-off False Discovery Rate below which to select DEGs
#' @param nom_pval the cut-off nominal P-value below which to select DEGs (as an alternative to FDR)
#' @param Nperms number of subsets created when downsampling at each level
#' @param y the column name in the SCE object for the return variable e.g. "diagnosis" - Case or disease. Default is the last variable in the design formula. y can be discrete (logistic regression) or continuous (linear regression)
#' @param region the column name in the SCE object for the study region. If there are multiple regions in the study (for example two brain regions). Pseudobulk values can be derived separately. Default is "single_region" which will not split by region.
#' @param control character specifying which control level for the differential expression analysis e.g. in a case/control/other study use "control" in the y column to compare against. NOTE only need to specify if more than two groups in y, leave as default value for two groups or continuous y. Default is NULL.
#' @param pval_adjust_method the adjustment method for the p-value in the differential expression analysis. Default is benjamini hochberg "BH". See  stats::p.adjust for available options
#' @param rmv_zero_count_genes whether genes with no count values in any cell should be removed. Default is TRUE

#' Saves all plots in the appropriate directory

power_plots <- function(SCE,
                        range_downsampled="placeholder",
                        output_path=getwd(),
                        inpath="placeholder",
                        sampled="individuals",
                        sampleID="donor_id",
                        celltypeID="cell_type",
                        coeff="male",
                        fdr=0.05,
                        nom_pval=0.05,
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
    DEGs_fdr <- subset(DEout$celltype_all_genes[[celltype_name]], adj_pval<fdr)$name
    DEGs_pval <- subset(DEout$celltype_all_genes[[celltype_name]], PValue<nom_pval)$name

    ## variables needed for the plots below (FDR cut-off)
    DEGs_full_fdr <- subset(DEout$celltype_all_genes[[celltype_name]], adj_pval<fdr)
    allgenes_fdr <- DEout$celltype_all_genes[[celltype_name]]
    totalDEGs_fdr <- length(DEGs_fdr)
    DEGs_full_fdr_upreg <- subset(DEout$celltype_all_genes[[celltype_name]], adj_pval<fdr&logFC > 0)
    DEGs_full_fdr_downreg <- subset(DEout$celltype_all_genes[[celltype_name]], adj_pval<fdr&logFC < 0)
    # look at quantiles_fdr of |logFC| and pick 25%, 50%, 75%, round to 1dp
    eff_size_fdr <- abs(subset(DEout$celltype_all_genes[[celltype_name]], adj_pval<fdr)$logFC)
    quantiles_fdr <- quantile(eff_size_fdr)
    eff_size_fdr2 <- subset(DEout$celltype_all_genes[[celltype_name]], adj_pval<fdr)$logFC
    upreg_fdr <- eff_size_fdr2[eff_size_fdr2 > 0]
    downreg_fdr <- eff_size_fdr2[eff_size_fdr2 < 0]
    quantiles_fdr_upreg <- quantile(upreg_fdr)
    quantiles_fdr_downreg <- quantile(downreg_fdr)
    # define quantiles_fdr
    Q1_fdr <- round(quantiles_fdr[[2]],1)
    Q2_fdr <- round(quantiles_fdr[[3]],1)
    Q3_fdr <- round(quantiles_fdr[[4]],1)
    Q1_fdr_upreg <- round(quantiles_fdr_upreg[[2]],1)
    Q2_fdr_upreg <- round(quantiles_fdr_upreg[[3]],1)
    Q3_fdr_upreg <- round(quantiles_fdr_upreg[[4]],1)
    Q1_fdr_downreg <- round(quantiles_fdr_downreg[[2]],1)
    Q2_fdr_downreg <- round(quantiles_fdr_downreg[[3]],1)
    Q3_fdr_downreg <- round(quantiles_fdr_downreg[[4]],1)
    # numbers of DEGs in each logFC range
    below_Q1_fdr <- length(eff_size_fdr[eff_size_fdr < Q1_fdr])
    below_Q1_Q2_fdr <- length(eff_size_fdr[eff_size_fdr >= Q1_fdr & eff_size_fdr < Q2_fdr])
    below_Q2_Q3_fdr <- length(eff_size_fdr[eff_size_fdr >= Q2_fdr & eff_size_fdr < Q3_fdr])
    above_Q3_fdr <- length(eff_size_fdr[eff_size_fdr >= Q3_fdr])
    below_Q1_fdr_upreg <- length(upreg_fdr[upreg_fdr < Q1_fdr_upreg])
    below_Q1_Q2_fdr_upreg <- length(upreg_fdr[upreg_fdr >= Q1_fdr_upreg & upreg_fdr < Q2_fdr_upreg])
    below_Q2_Q3_fdr_upreg <- length(upreg_fdr[upreg_fdr >= Q2_fdr_upreg & upreg_fdr < Q3_fdr_upreg])
    above_Q3_fdr_upreg <- length(upreg_fdr[upreg_fdr >= Q3_fdr_upreg])
    below_Q1_fdr_downreg <- length(downreg_fdr[downreg_fdr < Q1_fdr_downreg])
    below_Q1_Q2_fdr_downreg <- length(downreg_fdr[downreg_fdr >= Q1_fdr_downreg & downreg_fdr < Q2_fdr_downreg])
    below_Q2_Q3_fdr_downreg <- length(downreg_fdr[downreg_fdr >= Q2_fdr_downreg & downreg_fdr < Q3_fdr_downreg])
    above_Q3_fdr_downreg <- length(downreg_fdr[downreg_fdr >= Q3_fdr_downreg])
    # names of DEGs in each logFC range
    DEGs_below_Q1_fdr <- subset(DEGs_full_fdr,abs(logFC) < Q1_fdr)$name
    DEGs_between_Q1_Q2_fdr <- subset(DEGs_full_fdr,abs(logFC) >= Q1_fdr & abs(logFC) < Q2_fdr)$name
    DEGs_between_Q2_Q3_fdr <- subset(DEGs_full_fdr,abs(logFC) >= Q2_fdr & abs(logFC) < Q3_fdr)$name
    DEGs_above_Q3_fdr <- subset(DEGs_full_fdr,abs(logFC) >= Q3_fdr)$name
    nonDEGs_fdr <- subset(allgenes_fdr,!name%in%DEGs_fdr)
    DEGs_below_Q1_fdr_upreg <- subset(DEGs_full_fdr_upreg,logFC < Q1_fdr_upreg)$name
    DEGs_between_Q1_Q2_fdr_upreg <- subset(DEGs_full_fdr_upreg,logFC >= Q1_fdr_upreg & logFC < Q2_fdr_upreg)$name
    DEGs_between_Q2_Q3_fdr_upreg <- subset(DEGs_full_fdr_upreg,logFC >= Q2_fdr_upreg & logFC < Q3_fdr_upreg)$name
    DEGs_above_Q3_fdr_upreg <- subset(DEGs_full_fdr_upreg,logFC >= Q3_fdr_upreg)$name
    DEGs_below_Q1_fdr_downreg <- subset(DEGs_full_fdr_downreg,logFC < Q1_fdr_downreg)$name
    DEGs_between_Q1_Q2_fdr_downreg <- subset(DEGs_full_fdr_downreg,logFC >= Q1_fdr_downreg & logFC < Q2_fdr_downreg)$name
    DEGs_between_Q2_Q3_fdr_downreg <- subset(DEGs_full_fdr_downreg,logFC >= Q2_fdr_downreg & logFC < Q3_fdr_downreg)$name
    DEGs_above_Q3_fdr_downreg <- subset(DEGs_full_fdr_downreg,logFC >= Q3_fdr_downreg)$name
    # names of nonDEGs_fdr in each logFC range
    nonDEGs_belowQ1_fdr <- subset(nonDEGs_fdr,abs(logFC) < Q1_fdr)$name
    nonDEGs_Q1_to_Q2_fdr <- subset(nonDEGs_fdr,abs(logFC) >= Q1_fdr&abs(logFC) < Q2_fdr)$name
    nonDEGs_Q2_to_Q3_fdr <- subset(nonDEGs_fdr,abs(logFC) >= Q2_fdr&abs(logFC) < Q3_fdr)$name
    nonDEGs_aboveQ3_fdr <- subset(nonDEGs_fdr,abs(logFC) > Q3_fdr)$name
    ## variables needed for the plots below (pval cut-off)
    DEGs_full_pval <- subset(DEout$celltype_all_genes[[celltype_name]], PValue<nom_pval)
    allgenes_pval <- DEout$celltype_all_genes[[celltype_name]]
    totalDEGs_pval <- length(DEGs_pval)
    DEGs_full_pval_upreg <- subset(DEout$celltype_all_genes[[celltype_name]], PValue<nom_pval&logFC > 0)
    DEGs_full_pval_downreg <- subset(DEout$celltype_all_genes[[celltype_name]], PValue<nom_pval&logFC < 0)
    # look at quantiles_pval of |logFC| and pick 25%, 50%, 75%, round to 1dp
    eff_size_pval <- abs(subset(DEout$celltype_all_genes[[celltype_name]], PValue<nom_pval)$logFC)
    quantiles_pval <- quantile(eff_size_pval)
    eff_size_pval2 <- subset(DEout$celltype_all_genes[[celltype_name]], PValue<nom_pval)$logFC
    upreg_pval <- eff_size_pval2[eff_size_pval2 > 0]
    downreg_pval <- eff_size_pval2[eff_size_pval2 < 0]
    quantiles_pval_upreg <- quantile(upreg_pval)
    quantiles_pval_downreg <- quantile(downreg_pval)
    # define quantiles_pval
    Q1_pval <- round(quantiles_pval[[2]],1)
    Q2_pval <- round(quantiles_pval[[3]],1)
    Q3_pval <- round(quantiles_pval[[4]],1)
    Q1_pval_upreg <- round(quantiles_pval_upreg[[2]],1)
    Q2_pval_upreg <- round(quantiles_pval_upreg[[3]],1)
    Q3_pval_upreg <- round(quantiles_pval_upreg[[4]],1)
    Q1_pval_downreg <- round(quantiles_pval_downreg[[2]],1)
    Q2_pval_downreg <- round(quantiles_pval_downreg[[3]],1)
    Q3_pval_downreg <- round(quantiles_pval_downreg[[4]],1)
    # numbers of DEGs in each logFC range
    below_Q1_pval <- length(eff_size_pval[eff_size_pval < Q1_pval])
    below_Q1_Q2_pval <- length(eff_size_pval[eff_size_pval >= Q1_pval & eff_size_pval < Q2_pval])
    below_Q2_Q3_pval <- length(eff_size_pval[eff_size_pval >= Q2_pval & eff_size_pval < Q3_pval])
    above_Q3_pval <- length(eff_size_pval[eff_size_pval >= Q3_pval])
    below_Q1_pval_upreg <- length(upreg_pval[upreg_pval < Q1_pval_upreg])
    below_Q1_Q2_pval_upreg <- length(upreg_pval[upreg_pval >= Q1_pval_upreg & upreg_pval < Q2_pval_upreg])
    below_Q2_Q3_pval_upreg <- length(upreg_pval[upreg_pval >= Q2_pval_upreg & upreg_pval < Q3_pval_upreg])
    above_Q3_pval_upreg <- length(upreg_pval[upreg_pval >= Q3_pval_upreg])
    below_Q1_pval_downreg <- length(downreg_pval[downreg_pval < Q1_pval_downreg])
    below_Q1_Q2_pval_downreg <- length(downreg_pval[downreg_pval >= Q1_pval_downreg & downreg_pval < Q2_pval_downreg])
    below_Q2_Q3_pval_downreg <- length(downreg_pval[downreg_pval >= Q2_pval_downreg & downreg_pval < Q3_pval_downreg])
    above_Q3_pval_downreg <- length(downreg_pval[downreg_pval >= Q3_pval_downreg])
    # names of DEGs in each logFC range
    DEGs_below_Q1_pval <- subset(DEGs_full_pval,abs(logFC) < Q1_pval)$name
    DEGs_between_Q1_Q2_pval <- subset(DEGs_full_pval,abs(logFC) >= Q1_pval & abs(logFC) < Q2_pval)$name
    DEGs_between_Q2_Q3_pval <- subset(DEGs_full_pval,abs(logFC) >= Q2_pval & abs(logFC) < Q3_pval)$name
    DEGs_above_Q3_pval <- subset(DEGs_full_pval,abs(logFC) >= Q3_pval)$name
    nonDEGs_pval <- subset(allgenes_pval,!name%in%DEGs_pval)
    DEGs_below_Q1_pval_upreg <- subset(DEGs_full_pval_upreg,logFC < Q1_pval_upreg)$name
    DEGs_between_Q1_Q2_pval_upreg <- subset(DEGs_full_pval_upreg,logFC >= Q1_pval_upreg & logFC < Q2_pval_upreg)$name
    DEGs_between_Q2_Q3_pval_upreg <- subset(DEGs_full_pval_upreg,logFC >= Q2_pval_upreg & logFC < Q3_pval_upreg)$name
    DEGs_above_Q3_pval_upreg <- subset(DEGs_full_pval_upreg,logFC >= Q3_pval_upreg)$name
    DEGs_below_Q1_pval_downreg <- subset(DEGs_full_pval_downreg,logFC < Q1_pval_downreg)$name
    DEGs_between_Q1_Q2_pval_downreg <- subset(DEGs_full_pval_downreg,logFC >= Q1_pval_downreg & logFC < Q2_pval_downreg)$name
    DEGs_between_Q2_Q3_pval_downreg <- subset(DEGs_full_pval_downreg,logFC >= Q2_pval_downreg & logFC < Q3_pval_downreg)$name
    DEGs_above_Q3_pval_downreg <- subset(DEGs_full_pval_downreg,logFC >= Q3_pval_downreg)$name
    # names of nonDEGs_pval in each logFC range
    nonDEGs_belowQ1_pval <- subset(nonDEGs_pval,abs(logFC) < Q1_pval)$name
    nonDEGs_Q1_to_Q2_pval <- subset(nonDEGs_pval,abs(logFC) >= Q1_pval&abs(logFC) < Q2_pval)$name
    nonDEGs_Q2_to_Q3_pval <- subset(nonDEGs_pval,abs(logFC) >= Q2_pval&abs(logFC) < Q3_pval)$name
    nonDEGs_aboveQ3_pval <- subset(nonDEGs_pval,abs(logFC) > Q3_pval)$name

    if(sampled=="individuals"){
        
        #### plots using DEGs selected with an FDR cut-off
        # define path and folders
        path <- file.path(inpath,"DE_downsampling/")
        savepath <- file.path(output_path,"DE_downsampling/")
        dir.create(savepath,showWarnings=FALSE,recursive=TRUE)
        downsampled_folders <- paste0(paste(range_downsampled, sep=" ", collapse=NULL),"samples")
        ## plot % DEGs detected
        # load df with number of DEGs for each iteration/number of samples
        load(file.path(path, "DEGs_detected_fdr.RData"))
        # remove iteration column
        DEGs_detected_fdr <- subset(DEGs_detected_fdr,select=-c(Iteration))
        # convert to percentage
        DEGs_detected_fdr <- round(DEGs_detected_fdr/totalDEGs_fdr * 100, 2)
        # rename columns
        colnames(DEGs_detected_fdr) <- range_downsampled
        # reshape data
        DEGs_detected_new_fdr <- melt(DEGs_detected_fdr)
        ## create and save boxplots plot
        degs_detected_boxplot_fdr.plot<-
            ggplot(data=DEGs_detected_new_fdr,
                    aes(x = factor(variable), y = value))+
            geom_boxplot(fill="#69b3a2",color="black")+
            ggtitle("Percentage of total DEGs detected when down-sampling at a range of values (FDR)")+
            labs(y = "% DEGs detected", x = "Number of samples", fill="Samples") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9),legend.position="none")+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "DEGs_detected_boxplot_fdr.png"),degs_detected_boxplot_fdr.plot,width=25,height=20,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_detected_boxplot_fdr.pdf"),degs_detected_boxplot_fdr.plot,width=25,height=20,units="cm",bg="white")

        ## plot % of DEGs detected in each effect size range
        setwd(path)
        # loop through directories
        allDEGs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                DEGcounts <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], adj_pval<fdr)
                    # add numDEGs for each permutation
                    DEGcounts <- unlist(c(DEGcounts,
                                        length(subset(degs,name%in%DEGs_below_Q1_fdr)$name),
                                        length(subset(degs,name%in%DEGs_between_Q1_Q2_fdr)$name),
                                        length(subset(degs,name%in%DEGs_between_Q2_Q3_fdr)$name),
                                        length(subset(degs,name%in%DEGs_above_Q3_fdr)$name)))
                }
                allDEGs <- unlist(c(allDEGs, DEGcounts))
            }
        }
        # split by number of samples
        Iteration <- c()
        for(perm in 1:Nperms){
            Iteration <- c(Iteration, paste0(toString(perm),"_1"),
                                      paste0(toString(perm),"_2"),
                                      paste0(toString(perm),"_3"),
                                      paste0(toString(perm),"_4"))
        }
        # create dataframe
        allDEGs_effects <- data.frame(allDEGs[1:(4*Nperms)])
        colnames(allDEGs_effects) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            allDEGs_effects[[toString(range_downsampled[[value]])]] <- allDEGs[(4*(value-1)*Nperms+1):(4*value*Nperms)]
        }
        # convert to percentages
        for(j in 1:dim(allDEGs_effects)[[1]]){
            if(j%%4==1){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_fdr*100,2)
            }else if(j%%4==2){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_Q2_fdr*100,2)
            }else if(j%%4==3){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q2_Q3_fdr*100,2)
            }else{
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/above_Q3_fdr*100,2)
            }
        }
        # add and change iteration column above
        allDEGs_effects$`Iteration` <- Iteration
        allDEGs_effects$Iteration <- str_sub(allDEGs_effects$Iteration,-1,-1)
        # change name of iteration to be range of LogFCs
        allDEGs_effects$Iteration <- gsub("^1$",paste0("|logFC| < ",toString(Q1_fdr)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^2$",paste0(toString(Q1_fdr)," <= |logFC| < ",toString(Q2_fdr)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^3$",paste0(toString(Q2_fdr)," <= |logFC| < ",toString(Q3_fdr)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^4$",paste0("|logFC| >= ",toString(Q3_fdr)),allDEGs_effects$Iteration)
        # reshape
        allDEGs_effects_new <- melt(allDEGs_effects)
        # plot
        degs_effects_fdr.plot <- 
            ggplot(data=allDEGs_effects_new,
                aes(x = factor(variable), y = value, fill = factor(Iteration,levels=c(paste0("|logFC| < ",toString(Q1_fdr)),
                                                                                      paste0(toString(Q1_fdr)," <= |logFC| < ",toString(Q2_fdr)),
                                                                                      paste0(toString(Q2_fdr)," <= |logFC| < ",toString(Q3_fdr)),
                                                                                      paste0("|logFC| >= ",toString(Q3_fdr))),
                                                                             labels=c(paste0("|logFC| < ",toString(Q1_fdr)," (",below_Q1_fdr," DEGs)"),
                                                                                      paste0(toString(Q1_fdr)," <= |logFC| < ",toString(Q2_fdr)," (",below_Q1_Q2_fdr," DEGs)"),
                                                                                      paste0(toString(Q2_fdr)," <= |logFC| < ",toString(Q3_fdr)," (",below_Q2_Q3_fdr," DEGs)"),
                                                                                      paste0("|logFC| >= ",toString(Q3_fdr)," (",above_Q3_fdr," DEGs)")))))+
            geom_boxplot()+
            ggtitle("Percentage of total DEGs detected when down-sampling samples at a range of values, by effect sizes (FDR)")+
            labs(y = "% DEGs detected", x = "Number of samples", fill="Effect size range") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9))+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "DEGs_effects_fdr.png"),degs_effects_fdr.plot,width=35,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_effects_fdr.pdf"),degs_effects_fdr.plot,width=35,height=25,units="cm",bg="white")

        ## upreg plot
        setwd(path)
        # loop through directories
        allDEGs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                DEGcounts <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], adj_pval<fdr&logFC>0)
                    # add numDEGs for each permutation
                    DEGcounts <- unlist(c(DEGcounts,
                                        length(subset(degs,name%in%DEGs_below_Q1_fdr_upreg)$name),
                                        length(subset(degs,name%in%DEGs_between_Q1_Q2_fdr_upreg)$name),
                                        length(subset(degs,name%in%DEGs_between_Q2_Q3_fdr_upreg)$name),
                                        length(subset(degs,name%in%DEGs_above_Q3_fdr_upreg)$name)))
                }
                allDEGs <- unlist(c(allDEGs, DEGcounts))
            }
        }
        # split by number of samples
        Iteration <- c()
        for(perm in 1:Nperms){
            Iteration <- c(Iteration, paste0(toString(perm),"_1"),
                                    paste0(toString(perm),"_2"),
                                    paste0(toString(perm),"_3"),
                                    paste0(toString(perm),"_4"))
        }
        # create dataframe
        allDEGs_effects <- data.frame(allDEGs[1:(4*Nperms)])
        colnames(allDEGs_effects) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            allDEGs_effects[[toString(range_downsampled[[value]])]] <- allDEGs[(4*(value-1)*Nperms+1):(4*value*Nperms)]
        }
        # convert to percentages
        for(j in 1:dim(allDEGs_effects)[[1]]){
            if(j%%4==1){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_fdr_upreg*100,2)
            }else if(j%%4==2){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_Q2_fdr_upreg*100,2)
            }else if(j%%4==3){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q2_Q3_fdr_upreg*100,2)
            }else{
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/above_Q3_fdr_upreg*100,2)
            }
        }
        # add and change iteration column above
        allDEGs_effects$`Iteration` <- Iteration
        allDEGs_effects$Iteration <- str_sub(allDEGs_effects$Iteration,-1,-1)
        # change name of iteration to be range of LogFCs
        allDEGs_effects$Iteration <- gsub("^1$",paste0("0 < logFC < ",toString(Q1_fdr_upreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^2$",paste0(toString(Q1_fdr_upreg)," <= logFC < ",toString(Q2_fdr_upreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^3$",paste0(toString(Q2_fdr_upreg)," <= logFC < ",toString(Q3_fdr_upreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^4$",paste0("logFC >= ",toString(Q3_fdr_upreg)),allDEGs_effects$Iteration)
        # reshape
        allDEGs_effects_new <- melt(allDEGs_effects)
        # plot
        degs_effects_fdr_upreg.plot <- 
            ggplot(data=allDEGs_effects_new,
                aes(x = factor(variable), y = value, fill = factor(Iteration,levels=c(paste0("0 < logFC < ",toString(Q1_fdr_upreg)),
                                                                                        paste0(toString(Q1_fdr_upreg)," <= logFC < ",toString(Q2_fdr_upreg)),
                                                                                        paste0(toString(Q2_fdr_upreg)," <= logFC < ",toString(Q3_fdr_upreg)),
                                                                                        paste0("logFC >= ",toString(Q3_fdr_upreg))),
                                                                                labels=c(paste0("0 < logFC < ",toString(Q1_fdr_upreg)," (",below_Q1_fdr_upreg," DEGs)"),
                                                                                        paste0(toString(Q1_fdr_upreg)," <= logFC < ",toString(Q2_fdr_upreg)," (",below_Q1_Q2_fdr_upreg," DEGs)"),
                                                                                        paste0(toString(Q2_fdr_upreg)," <= logFC < ",toString(Q3_fdr_upreg)," (",below_Q2_Q3_fdr_upreg," DEGs)"),
                                                                                        paste0("logFC >= ",toString(Q3_fdr_upreg)," (",above_Q3_fdr_upreg," DEGs)")))))+
            geom_boxplot()+
            ggtitle("Percentage of total (up-reg) DEGs detected when down-sampling at a range of values, by effect sizes (FDR)")+
            labs(y = "% DEGs detected", x = "Number of samples", fill="Effect size range") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9))+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "DEGs_upreg_effects_fdr.png"),degs_effects_fdr_upreg.plot,width=35,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_upreg_effects_fdr.pdf"),degs_effects_fdr_upreg.plot,width=35,height=25,units="cm",bg="white")
        
        ## downreg plot
        setwd(path)
        # loop through directories
        allDEGs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                DEGcounts <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], adj_pval<fdr&logFC<0)
                    # add numDEGs for each permutation
                    DEGcounts <- unlist(c(DEGcounts,
                                        length(subset(degs,name%in%DEGs_below_Q1_fdr_downreg)$name),
                                        length(subset(degs,name%in%DEGs_between_Q1_Q2_fdr_downreg)$name),
                                        length(subset(degs,name%in%DEGs_between_Q2_Q3_fdr_downreg)$name),
                                        length(subset(degs,name%in%DEGs_above_Q3_fdr_downreg)$name)))
                }
                allDEGs <- unlist(c(allDEGs, DEGcounts))
            }
        }
        # split by number of samples
        Iteration <- c()
        for(perm in 1:Nperms){
            Iteration <- c(Iteration, paste0(toString(perm),"_1"),
                                    paste0(toString(perm),"_2"),
                                    paste0(toString(perm),"_3"),
                                    paste0(toString(perm),"_4"))
        }
        # create dataframe
        allDEGs_effects <- data.frame(allDEGs[1:(4*Nperms)])
        colnames(allDEGs_effects) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            allDEGs_effects[[toString(range_downsampled[[value]])]] <- allDEGs[(4*(value-1)*Nperms+1):(4*value*Nperms)]
        }
        # convert to percentages
        for(j in 1:dim(allDEGs_effects)[[1]]){
            if(j%%4==1){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_fdr_downreg*100,2)
            }else if(j%%4==2){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_Q2_fdr_downreg*100,2)
            }else if(j%%4==3){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q2_Q3_fdr_downreg*100,2)
            }else{
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/above_Q3_fdr_downreg*100,2)
            }
        }
        # add and change iteration column above
        allDEGs_effects$`Iteration` <- Iteration
        allDEGs_effects$Iteration <- str_sub(allDEGs_effects$Iteration,-1,-1)
        # change name of iteration to be range of LogFCs
        allDEGs_effects$Iteration <- gsub("^1$",paste0("logFC < ",toString(Q1_fdr_downreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^2$",paste0(toString(Q1_fdr_downreg)," <= logFC < ",toString(Q2_fdr_downreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^3$",paste0(toString(Q2_fdr_downreg)," <= logFC < ",toString(Q3_fdr_downreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^4$",paste0("0 > logFC >= ",toString(Q3_fdr_downreg)),allDEGs_effects$Iteration)
        # reshape
        allDEGs_effects_new <- melt(allDEGs_effects)
        # plot
        degs_effects_fdr_downreg.plot <- 
            ggplot(data=allDEGs_effects_new,
                aes(x = factor(variable), y = value, fill = factor(Iteration,levels=c(paste0("logFC < ",toString(Q1_fdr_downreg)),
                                                                                        paste0(toString(Q1_fdr_downreg)," <= logFC < ",toString(Q2_fdr_downreg)),
                                                                                        paste0(toString(Q2_fdr_downreg)," <= logFC < ",toString(Q3_fdr_downreg)),
                                                                                        paste0("0 > logFC >= ",toString(Q3_fdr_downreg))),
                                                                                labels=c(paste0("logFC < ",toString(Q1_fdr_downreg)," (",below_Q1_fdr_downreg," DEGs)"),
                                                                                        paste0(toString(Q1_fdr_downreg)," <= logFC < ",toString(Q2_fdr_downreg)," (",below_Q1_Q2_fdr_downreg," DEGs)"),
                                                                                        paste0(toString(Q2_fdr_downreg)," <= logFC < ",toString(Q3_fdr_downreg)," (",below_Q2_Q3_fdr_downreg," DEGs)"),
                                                                                        paste0("0 > logFC >= ",toString(Q3_fdr_downreg)," (",above_Q3_fdr_downreg," DEGs)")))))+
            geom_boxplot()+
            ggtitle("Percentage of total (down-reg) DEGs detected when down-sampling at a range of values, by effect sizes (FDR)")+
            labs(y = "% DEGs detected", x = "Number of samples", fill="Effect size range") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9))+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "DEGs_downreg_effects_fdr.png"),degs_effects_fdr_downreg.plot,width=35,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_downreg_effects_fdr.pdf"),degs_effects_fdr_downreg.plot,width=35,height=25,units="cm",bg="white")

        ## plot False Positive Rate
        setwd(path)
        # loop through directories
        allFPRs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                FPRs <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], adj_pval<fdr)$name
                    allgenes_tmp <- get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]]
                    non_degs <- subset(allgenes_tmp,!name%in%degs)
                    TNs <- length(subset(non_degs,name%in%nonDEGs_fdr$name)$name)
                    # get FPR
                    numFPs <- length(degs[degs %in% nonDEGs_fdr$name])
                    FPR <- numFPs/(numFPs+TNs)
                    # add numFPs for each permutation
                    FPRs <- unlist(c(FPRs, FPR))
                }
                allFPRs <- unlist(c(allFPRs, FPRs))
            }
        }
        # create dataframe
        all_FPRs <- data.frame(allFPRs[1:Nperms])
        colnames(all_FPRs) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            all_FPRs[[toString(range_downsampled[[value]])]] <- allFPRs[((value-1)*Nperms+1):(value*Nperms)]
        }
        # reshape data
        all_FPRs_new <- melt(all_FPRs)
        # create and save boxplots
        FPs_detected_boxplot_fdr.plot<-
            ggplot(data=all_FPRs_new,
                    aes(x = factor(variable), y = value))+
            geom_boxplot(fill="#f37735",color="black")+
            ggtitle("False Positive Rate (FPR) when down-sampling at a range of values (FDR)")+
            labs(y = "FPR", x = "Number of samples", fill="Samples") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9),legend.position="none")+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "FPs_detected_boxplot_fdr.png"),FPs_detected_boxplot_fdr.plot,width=25,height=20,units="cm",bg="white")
        ggsave(file.path(savepath, "FPs_detected_boxplot_fdr.pdf"),FPs_detected_boxplot_fdr.plot,width=25,height=20,units="cm",bg="white")

        ## plot False Discovery Rate
        setwd(path)
        # loop through directories
        allFDRs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                FDRs <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], adj_pval<fdr)$name
                    # get true positives
                    numTPs <- length(degs[degs %in% DEGs_fdr])
                    # get FDR
                    numFPs <- length(degs[degs %in% nonDEGs_fdr$name])
                    FDR <- numFPs/(numFPs+numTPs)
                    # add FDR for each permutation
                    FDRs <- unlist(c(FDRs, FDR))
                }
                allFDRs <- unlist(c(allFDRs, FDRs))
            }
        }
        # create dataframe
        all_FDRs <- data.frame(allFDRs[1:Nperms])
        colnames(all_FDRs) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            all_FDRs[[toString(range_downsampled[[value]])]] <- allFDRs[((value-1)*Nperms+1):(value*Nperms)]
        }
        # reshape data
        all_FDRs_new <- melt(all_FDRs)
        # create and save boxplots
        FDRs_boxplot_fdr.plot<-
            ggplot(data=all_FDRs_new,
                    aes(x = factor(variable), y = value))+
            geom_boxplot(fill="#f37735",color="black")+
            ggtitle("False Discovery Rate (FDR) when down-sampling at a range of values (FDR)")+
            labs(y = "FDR", x = "Number of samples", fill="Samples") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9),legend.position="none")+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "FDRs_boxplot_fdr.png"),FDRs_boxplot_fdr.plot,width=25,height=20,units="cm",bg="white")
        ggsave(file.path(savepath, "FDRs_boxplot_fdr.pdf"),FDRs_boxplot_fdr.plot,width=25,height=20,units="cm",bg="white")        

        ## plot False Discovery Rate of DEGs in each effect size range
        setwd(path)
        # loop through directories
        allFDRs_effs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                FDRs <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], adj_pval<fdr)
                    # get numTPs for each range
                    TP_belowQ1 <- length(subset(degs,name%in%DEGs_below_Q1_fdr)$name)
                    TP_Q1_to_Q2 <- length(subset(degs,name%in%DEGs_between_Q1_Q2_fdr)$name)
                    TP_Q2_to_Q3 <- length(subset(degs,name%in%DEGs_between_Q2_Q3_fdr)$name)
                    TP_aboveQ3 <- length(subset(degs,name%in%DEGs_above_Q3_fdr)$name)
                    # get FPs for each range
                    FP_belowQ1 <- length(subset(degs,name%in%nonDEGs_belowQ1_fdr)$name)
                    FP_Q1_to_Q2 <- length(subset(degs,name%in%nonDEGs_Q1_to_Q2_fdr)$name)
                    FP_Q2_to_Q3 <- length(subset(degs,name%in%nonDEGs_Q2_to_Q3_fdr)$name)
                    FP_aboveQ3 <- length(subset(degs,name%in%nonDEGs_aboveQ3_fdr)$name)
                    # add numDEGs for each permutation
                    FDRs <- unlist(c(FDRs,
                                    FP_belowQ1/(FP_belowQ1+TP_belowQ1),
                                    FP_Q1_to_Q2/(FP_Q1_to_Q2+TP_Q1_to_Q2),
                                    FP_Q2_to_Q3/(FP_Q2_to_Q3+TP_Q2_to_Q3),
                                    FP_aboveQ3/(FP_aboveQ3+TP_aboveQ3)))
                }
                allFDRs_effs <- unlist(c(allFDRs_effs, FDRs))
            }
        }
        # split by number of samples
        Iteration <- c()
        for(perm in 1:Nperms){
            Iteration <- c(Iteration, paste0(toString(perm),"_1"),
                                      paste0(toString(perm),"_2"),
                                      paste0(toString(perm),"_3"),
                                      paste0(toString(perm),"_4"))
        }
        # create dataframe
        all_FDR_effs <- data.frame(allFDRs_effs[1:(4*Nperms)])
        colnames(all_FDR_effs) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            all_FDR_effs[[toString(range_downsampled[[value]])]] <- allFDRs_effs[(4*(value-1)*Nperms+1):(4*value*Nperms)]
        }
        ## create boxplots for various ranges of effect sizes
        # add and change iteration column above
        all_FDR_effs$`Iteration` <- Iteration
        all_FDR_effs$Iteration <- str_sub(all_FDR_effs$Iteration,-1,-1)
        # change name of iteration to be range of LogFCs
        all_FDR_effs$Iteration <- gsub("^1$",paste0("|logFC| < ",toString(Q1_fdr)),all_FDR_effs$Iteration)
        all_FDR_effs$Iteration <- gsub("^2$",paste0(toString(Q1_fdr)," <= |logFC| < ",toString(Q2_fdr)),all_FDR_effs$Iteration)
        all_FDR_effs$Iteration <- gsub("^3$",paste0(toString(Q2_fdr)," <= |logFC| < ",toString(Q3_fdr)),all_FDR_effs$Iteration)
        all_FDR_effs$Iteration <- gsub("^4$",paste0("|logFC| >= ",toString(Q3_fdr)),all_FDR_effs$Iteration)
        # reshape
        all_FDR_effs_new <- melt(all_FDR_effs)
        # plot
        FPs_effects_fdr.plot <- 
            ggplot(data=all_FDR_effs_new,
                aes(x = factor(variable), y = value, fill = factor(Iteration,levels=c(paste0("|logFC| < ",toString(Q1_fdr)),
                                                                                      paste0(toString(Q1_fdr)," <= |logFC| < ",toString(Q2_fdr)),
                                                                                      paste0(toString(Q2_fdr)," <= |logFC| < ",toString(Q3_fdr)),
                                                                                      paste0("|logFC| >= ",toString(Q3_fdr))),
                                                                             labels=c(paste0("|logFC| < ",toString(Q1_fdr)," (",below_Q1_fdr," DEGs)"),
                                                                                      paste0(toString(Q1_fdr)," <= |logFC| < ",toString(Q2_fdr)," (",below_Q1_Q2_fdr," DEGs)"),
                                                                                      paste0(toString(Q2_fdr)," <= |logFC| < ",toString(Q3_fdr)," (",below_Q2_Q3_fdr," DEGs)"),
                                                                                      paste0("|logFC| >= ",toString(Q3_fdr)," (",above_Q3_fdr," DEGs)")))))+
            geom_boxplot()+
            ggtitle("False Discovery Rate (FDR) when down-sampling samples at a range of values, by effect sizes (FDR)")+
            labs(y = "FDR", x = "Number of samples", fill="Effect size range") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9))+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "FPs_effects_fdr.png"),FPs_effects_fdr.plot,width=30,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "FPs_effects_fdr.pdf"),FPs_effects_fdr.plot,width=30,height=25,units="cm",bg="white")

        #### plots using DEGs selected with a nominal p-value cut-off
        ## plot % DEGs detected
        # load df with number of DEGs for each iteration/number of samples
        load(file.path(path, "DEGs_detected_pval.RData"))
        # remove iteration column
        DEGs_detected_pval <- subset(DEGs_detected_pval,select=-c(Iteration))
        # convert to percentage
        DEGs_detected_pval <- round(DEGs_detected_pval/totalDEGs_pval * 100, 2)
        # rename columns
        colnames(DEGs_detected_pval) <- range_downsampled
        # reshape data
        DEGs_detected_new_pval <- melt(DEGs_detected_pval)
        ## create and save boxplots plot
        degs_detected_boxplot_pval.plot<-
            ggplot(data=DEGs_detected_new_pval,
                    aes(x = factor(variable), y = value))+
            geom_boxplot(fill="#69b3a2",color="black")+
            ggtitle("Percentage of total DEGs detected when down-sampling at a range of values (p-val)")+
            labs(y = "% DEGs detected", x = "Number of samples", fill="Samples") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9),legend.position="none")+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "DEGs_detected_boxplot_pval.png"),degs_detected_boxplot_pval.plot,width=25,height=20,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_detected_boxplot_pval.pdf"),degs_detected_boxplot_pval.plot,width=25,height=20,units="cm",bg="white")

        ## plot % of DEGs detected in each effect size range
        setwd(path)
        # loop through directories
        allDEGs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                DEGcounts <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], PValue<nom_pval)
                    # add numDEGs for each permutation
                    DEGcounts <- unlist(c(DEGcounts,
                                        length(subset(degs,name%in%DEGs_below_Q1_pval)$name),
                                        length(subset(degs,name%in%DEGs_between_Q1_Q2_pval)$name),
                                        length(subset(degs,name%in%DEGs_between_Q2_Q3_pval)$name),
                                        length(subset(degs,name%in%DEGs_above_Q3_pval)$name)))
                }
                allDEGs <- unlist(c(allDEGs, DEGcounts))
            }
        }
        # split by number of samples
        Iteration <- c()
        for(perm in 1:Nperms){
            Iteration <- c(Iteration, paste0(toString(perm),"_1"),
                                      paste0(toString(perm),"_2"),
                                      paste0(toString(perm),"_3"),
                                      paste0(toString(perm),"_4"))
        }
        # create dataframe
        allDEGs_effects <- data.frame(allDEGs[1:(4*Nperms)])
        colnames(allDEGs_effects) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            allDEGs_effects[[toString(range_downsampled[[value]])]] <- allDEGs[(4*(value-1)*Nperms+1):(4*value*Nperms)]
        }
        # convert to percentages
        for(j in 1:dim(allDEGs_effects)[[1]]){
            if(j%%4==1){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_pval*100,2)
            }else if(j%%4==2){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_Q2_pval*100,2)
            }else if(j%%4==3){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q2_Q3_pval*100,2)
            }else{
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/above_Q3_pval*100,2)
            }
        }
        # add and change iteration column above
        allDEGs_effects$`Iteration` <- Iteration
        allDEGs_effects$Iteration <- str_sub(allDEGs_effects$Iteration,-1,-1)
        # change name of iteration to be range of LogFCs
        allDEGs_effects$Iteration <- gsub("^1$",paste0("|logFC| < ",toString(Q1_pval)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^2$",paste0(toString(Q1_pval)," <= |logFC| < ",toString(Q2_pval)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^3$",paste0(toString(Q2_pval)," <= |logFC| < ",toString(Q3_pval)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^4$",paste0("|logFC| >= ",toString(Q3_pval)),allDEGs_effects$Iteration)
        # reshape
        allDEGs_effects_new <- melt(allDEGs_effects)
        # plot
        degs_effects_pval.plot <- 
            ggplot(data=allDEGs_effects_new,
                aes(x = factor(variable), y = value, fill = factor(Iteration,levels=c(paste0("|logFC| < ",toString(Q1_pval)),
                                                                                      paste0(toString(Q1_pval)," <= |logFC| < ",toString(Q2_pval)),
                                                                                      paste0(toString(Q2_pval)," <= |logFC| < ",toString(Q3_pval)),
                                                                                      paste0("|logFC| >= ",toString(Q3_pval))),
                                                                             labels=c(paste0("|logFC| < ",toString(Q1_pval)," (",below_Q1_pval," DEGs)"),
                                                                                      paste0(toString(Q1_pval)," <= |logFC| < ",toString(Q2_pval)," (",below_Q1_Q2_pval," DEGs)"),
                                                                                      paste0(toString(Q2_pval)," <= |logFC| < ",toString(Q3_pval)," (",below_Q2_Q3_pval," DEGs)"),
                                                                                      paste0("|logFC| >= ",toString(Q3_pval)," (",above_Q3_pval," DEGs)")))))+
            geom_boxplot()+
            ggtitle("Percentage of total DEGs detected when down-sampling samples at a range of values, by effect sizes (p-val)")+
            labs(y = "% DEGs detected", x = "Number of samples", fill="Effect size range") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9))+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "DEGs_effects_pval.png"),degs_effects_pval.plot,width=35,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_effects_pval.pdf"),degs_effects_pval.plot,width=35,height=25,units="cm",bg="white")

        ## upreg plot
        setwd(path)
        # loop through directories
        allDEGs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                DEGcounts <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], PValue<nom_pval&logFC>0)
                    # add numDEGs for each permutation
                    DEGcounts <- unlist(c(DEGcounts,
                                        length(subset(degs,name%in%DEGs_below_Q1_pval_upreg)$name),
                                        length(subset(degs,name%in%DEGs_between_Q1_Q2_pval_upreg)$name),
                                        length(subset(degs,name%in%DEGs_between_Q2_Q3_pval_upreg)$name),
                                        length(subset(degs,name%in%DEGs_above_Q3_pval_upreg)$name)))
                }
                allDEGs <- unlist(c(allDEGs, DEGcounts))
            }
        }
        # split by number of samples
        Iteration <- c()
        for(perm in 1:Nperms){
            Iteration <- c(Iteration, paste0(toString(perm),"_1"),
                                    paste0(toString(perm),"_2"),
                                    paste0(toString(perm),"_3"),
                                    paste0(toString(perm),"_4"))
        }
        # create dataframe
        allDEGs_effects <- data.frame(allDEGs[1:(4*Nperms)])
        colnames(allDEGs_effects) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            allDEGs_effects[[toString(range_downsampled[[value]])]] <- allDEGs[(4*(value-1)*Nperms+1):(4*value*Nperms)]
        }
        # convert to percentages
        for(j in 1:dim(allDEGs_effects)[[1]]){
            if(j%%4==1){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_pval_upreg*100,2)
            }else if(j%%4==2){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_Q2_pval_upreg*100,2)
            }else if(j%%4==3){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q2_Q3_pval_upreg*100,2)
            }else{
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/above_Q3_pval_upreg*100,2)
            }
        }
        # add and change iteration column above
        allDEGs_effects$`Iteration` <- Iteration
        allDEGs_effects$Iteration <- str_sub(allDEGs_effects$Iteration,-1,-1)
        # change name of iteration to be range of LogFCs
        allDEGs_effects$Iteration <- gsub("^1$",paste0("0 < logFC < ",toString(Q1_pval_upreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^2$",paste0(toString(Q1_pval_upreg)," <= logFC < ",toString(Q2_pval_upreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^3$",paste0(toString(Q2_pval_upreg)," <= logFC < ",toString(Q3_pval_upreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^4$",paste0("logFC >= ",toString(Q3_pval_upreg)),allDEGs_effects$Iteration)
        # reshape
        allDEGs_effects_new <- melt(allDEGs_effects)
        # plot
        degs_effects_pval_upreg.plot <- 
            ggplot(data=allDEGs_effects_new,
                aes(x = factor(variable), y = value, fill = factor(Iteration,levels=c(paste0("0 < logFC < ",toString(Q1_pval_upreg)),
                                                                                        paste0(toString(Q1_pval_upreg)," <= logFC < ",toString(Q2_pval_upreg)),
                                                                                        paste0(toString(Q2_pval_upreg)," <= logFC < ",toString(Q3_pval_upreg)),
                                                                                        paste0("logFC >= ",toString(Q3_pval_upreg))),
                                                                                labels=c(paste0("0 < logFC < ",toString(Q1_pval_upreg)," (",below_Q1_pval_upreg," DEGs)"),
                                                                                        paste0(toString(Q1_pval_upreg)," <= logFC < ",toString(Q2_pval_upreg)," (",below_Q1_Q2_pval_upreg," DEGs)"),
                                                                                        paste0(toString(Q2_pval_upreg)," <= logFC < ",toString(Q3_pval_upreg)," (",below_Q2_Q3_pval_upreg," DEGs)"),
                                                                                        paste0("logFC >= ",toString(Q3_pval_upreg)," (",above_Q3_pval_upreg," DEGs)")))))+
            geom_boxplot()+
            ggtitle("Percentage of total (up-reg) DEGs detected when down-sampling at a range of values, by effect sizes (p-val)")+
            labs(y = "% DEGs detected", x = "Number of samples", fill="Effect size range") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9))+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "DEGs_upreg_effects_pval.png"),degs_effects_pval_upreg.plot,width=35,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_upreg_effects_pval.pdf"),degs_effects_pval_upreg.plot,width=35,height=25,units="cm",bg="white")
        
        ## downreg plot
        setwd(path)
        # loop through directories
        allDEGs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                DEGcounts <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], PValue<nom_pval&logFC<0)
                    # add numDEGs for each permutation
                    DEGcounts <- unlist(c(DEGcounts,
                                        length(subset(degs,name%in%DEGs_below_Q1_pval_downreg)$name),
                                        length(subset(degs,name%in%DEGs_between_Q1_Q2_pval_downreg)$name),
                                        length(subset(degs,name%in%DEGs_between_Q2_Q3_pval_downreg)$name),
                                        length(subset(degs,name%in%DEGs_above_Q3_pval_downreg)$name)))
                }
                allDEGs <- unlist(c(allDEGs, DEGcounts))
            }
        }
        # split by number of samples
        Iteration <- c()
        for(perm in 1:Nperms){
            Iteration <- c(Iteration, paste0(toString(perm),"_1"),
                                    paste0(toString(perm),"_2"),
                                    paste0(toString(perm),"_3"),
                                    paste0(toString(perm),"_4"))
        }
        # create dataframe
        allDEGs_effects <- data.frame(allDEGs[1:(4*Nperms)])
        colnames(allDEGs_effects) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            allDEGs_effects[[toString(range_downsampled[[value]])]] <- allDEGs[(4*(value-1)*Nperms+1):(4*value*Nperms)]
        }
        # convert to percentages
        for(j in 1:dim(allDEGs_effects)[[1]]){
            if(j%%4==1){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_pval_downreg*100,2)
            }else if(j%%4==2){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_Q2_pval_downreg*100,2)
            }else if(j%%4==3){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q2_Q3_pval_downreg*100,2)
            }else{
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/above_Q3_pval_downreg*100,2)
            }
        }
        # add and change iteration column above
        allDEGs_effects$`Iteration` <- Iteration
        allDEGs_effects$Iteration <- str_sub(allDEGs_effects$Iteration,-1,-1)
        # change name of iteration to be range of LogFCs
        allDEGs_effects$Iteration <- gsub("^1$",paste0("logFC < ",toString(Q1_pval_downreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^2$",paste0(toString(Q1_pval_downreg)," <= logFC < ",toString(Q2_pval_downreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^3$",paste0(toString(Q2_pval_downreg)," <= logFC < ",toString(Q3_pval_downreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^4$",paste0("0 > logFC >= ",toString(Q3_pval_downreg)),allDEGs_effects$Iteration)
        # reshape
        allDEGs_effects_new <- melt(allDEGs_effects)
        # plot
        degs_effects_pval_downreg.plot <- 
            ggplot(data=allDEGs_effects_new,
                aes(x = factor(variable), y = value, fill = factor(Iteration,levels=c(paste0("logFC < ",toString(Q1_pval_downreg)),
                                                                                        paste0(toString(Q1_pval_downreg)," <= logFC < ",toString(Q2_pval_downreg)),
                                                                                        paste0(toString(Q2_pval_downreg)," <= logFC < ",toString(Q3_pval_downreg)),
                                                                                        paste0("0 > logFC >= ",toString(Q3_pval_downreg))),
                                                                                labels=c(paste0("logFC < ",toString(Q1_pval_downreg)," (",below_Q1_pval_downreg," DEGs)"),
                                                                                        paste0(toString(Q1_pval_downreg)," <= logFC < ",toString(Q2_pval_downreg)," (",below_Q1_Q2_pval_downreg," DEGs)"),
                                                                                        paste0(toString(Q2_pval_downreg)," <= logFC < ",toString(Q3_pval_downreg)," (",below_Q2_Q3_pval_downreg," DEGs)"),
                                                                                        paste0("0 > logFC >= ",toString(Q3_pval_downreg)," (",above_Q3_pval_downreg," DEGs)")))))+
            geom_boxplot()+
            ggtitle("Percentage of total (down-reg) DEGs detected when down-sampling at a range of values, by effect sizes (p-val)")+
            labs(y = "% DEGs detected", x = "Number of samples", fill="Effect size range") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9))+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "DEGs_downreg_effects_pval.png"),degs_effects_pval_downreg.plot,width=35,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_downreg_effects_pval.pdf"),degs_effects_pval_downreg.plot,width=35,height=25,units="cm",bg="white")

        ## plot False Positive Rate
        setwd(path)
        # loop through directories
        allFPRs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                FPRs <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], PValue<nom_pval)$name
                    allgenes_tmp <- get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]]
                    non_degs <- subset(allgenes_tmp,!name%in%degs)
                    TNs <- length(subset(non_degs,name%in%nonDEGs_pval$name)$name)
                    # get FPR
                    numFPs <- length(degs[degs %in% nonDEGs_pval$name])
                    FPR <- numFPs/(numFPs+TNs)
                    # add numFPs for each permutation
                    FPRs <- unlist(c(FPRs, FPR))
                }
                allFPRs <- unlist(c(allFPRs, FPRs))
            }
        }
        # create dataframe
        all_FPRs <- data.frame(allFPRs[1:Nperms])
        colnames(all_FPRs) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            all_FPRs[[toString(range_downsampled[[value]])]] <- allFPRs[((value-1)*Nperms+1):(value*Nperms)]
        }
        # reshape data
        all_FPRs_new <- melt(all_FPRs)
        # create and save boxplots
        FPs_detected_boxplot_pval.plot<-
            ggplot(data=all_FPRs_new,
                    aes(x = factor(variable), y = value))+
            geom_boxplot(fill="#f37735",color="black")+
            ggtitle("False Positive Rate (FPR) when down-sampling at a range of values (p-val)")+
            labs(y = "FPR", x = "Number of samples", fill="Samples") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9),legend.position="none")+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "FPs_detected_boxplot_pval.png"),FPs_detected_boxplot_pval.plot,width=25,height=20,units="cm",bg="white")
        ggsave(file.path(savepath, "FPs_detected_boxplot_pval.pdf"),FPs_detected_boxplot_pval.plot,width=25,height=20,units="cm",bg="white")

        ## plot False Discovery Rate
        setwd(path)
        # loop through directories
        allFDRs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                FDRs <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], PValue<nom_pval)$name
                    # get true positives
                    numTPs <- length(degs[degs %in% DEGs_pval])
                    # get FDR
                    numFPs <- length(degs[degs %in% nonDEGs_pval$name])
                    FDR <- numFPs/(numFPs+numTPs)
                    # add FDR for each permutation
                    FDRs <- unlist(c(FDRs, FDR))
                }
                allFDRs <- unlist(c(allFDRs, FDRs))
            }
        }
        # create dataframe
        all_FDRs <- data.frame(allFDRs[1:Nperms])
        colnames(all_FDRs) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            all_FDRs[[toString(range_downsampled[[value]])]] <- allFDRs[((value-1)*Nperms+1):(value*Nperms)]
        }
        # reshape data
        all_FDRs_new <- melt(all_FDRs)
        # create and save boxplots
        FDRs_boxplot_pval.plot<-
            ggplot(data=all_FDRs_new,
                    aes(x = factor(variable), y = value))+
            geom_boxplot(fill="#f37735",color="black")+
            ggtitle("False Discovery Rate (FDR) when down-sampling at a range of values (p-val)")+
            labs(y = "FDR", x = "Number of samples", fill="Samples") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9),legend.position="none")+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "FDRs_boxplot_pval.png"),FDRs_boxplot_pval.plot,width=25,height=20,units="cm",bg="white")
        ggsave(file.path(savepath, "FDRs_boxplot_pval.pdf"),FDRs_boxplot_pval.plot,width=25,height=20,units="cm",bg="white")        

        ## plot False Discovery Rate of DEGs in each effect size range
        setwd(path)
        # loop through directories
        allFDRs_effs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                FDRs <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], PValue<nom_pval)
                    # get numTPs for each range
                    TP_belowQ1 <- length(subset(degs,name%in%DEGs_below_Q1_pval)$name)
                    TP_Q1_to_Q2 <- length(subset(degs,name%in%DEGs_between_Q1_Q2_pval)$name)
                    TP_Q2_to_Q3 <- length(subset(degs,name%in%DEGs_between_Q2_Q3_pval)$name)
                    TP_aboveQ3 <- length(subset(degs,name%in%DEGs_above_Q3_pval)$name)
                    # get FPs for each range
                    FP_belowQ1 <- length(subset(degs,name%in%nonDEGs_belowQ1_pval)$name)
                    FP_Q1_to_Q2 <- length(subset(degs,name%in%nonDEGs_Q1_to_Q2_pval)$name)
                    FP_Q2_to_Q3 <- length(subset(degs,name%in%nonDEGs_Q2_to_Q3_pval)$name)
                    FP_aboveQ3 <- length(subset(degs,name%in%nonDEGs_aboveQ3_pval)$name)
                    # add numDEGs for each permutation
                    FDRs <- unlist(c(FDRs,
                                    FP_belowQ1/(FP_belowQ1+TP_belowQ1),
                                    FP_Q1_to_Q2/(FP_Q1_to_Q2+TP_Q1_to_Q2),
                                    FP_Q2_to_Q3/(FP_Q2_to_Q3+TP_Q2_to_Q3),
                                    FP_aboveQ3/(FP_aboveQ3+TP_aboveQ3)))
                }
                allFDRs_effs <- unlist(c(allFDRs_effs, FDRs))
            }
        }
        # split by number of samples
        Iteration <- c()
        for(perm in 1:Nperms){
            Iteration <- c(Iteration, paste0(toString(perm),"_1"),
                                      paste0(toString(perm),"_2"),
                                      paste0(toString(perm),"_3"),
                                      paste0(toString(perm),"_4"))
        }
        # create dataframe
        all_FDR_effs <- data.frame(allFDRs_effs[1:(4*Nperms)])
        colnames(all_FDR_effs) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            all_FDR_effs[[toString(range_downsampled[[value]])]] <- allFDRs_effs[(4*(value-1)*Nperms+1):(4*value*Nperms)]
        }
        ## create boxplots for various ranges of effect sizes
        # add and change iteration column above
        all_FDR_effs$`Iteration` <- Iteration
        all_FDR_effs$Iteration <- str_sub(all_FDR_effs$Iteration,-1,-1)
        # change name of iteration to be range of LogFCs
        all_FDR_effs$Iteration <- gsub("^1$",paste0("|logFC| < ",toString(Q1_pval)),all_FDR_effs$Iteration)
        all_FDR_effs$Iteration <- gsub("^2$",paste0(toString(Q1_pval)," <= |logFC| < ",toString(Q2_pval)),all_FDR_effs$Iteration)
        all_FDR_effs$Iteration <- gsub("^3$",paste0(toString(Q2_pval)," <= |logFC| < ",toString(Q3_pval)),all_FDR_effs$Iteration)
        all_FDR_effs$Iteration <- gsub("^4$",paste0("|logFC| >= ",toString(Q3_pval)),all_FDR_effs$Iteration)
        # reshape
        all_FDR_effs_new <- melt(all_FDR_effs)
        # plot
        FPs_effects_pval.plot <- 
            ggplot(data=all_FDR_effs_new,
                aes(x = factor(variable), y = value, fill = factor(Iteration,levels=c(paste0("|logFC| < ",toString(Q1_pval)),
                                                                                      paste0(toString(Q1_pval)," <= |logFC| < ",toString(Q2_pval)),
                                                                                      paste0(toString(Q2_pval)," <= |logFC| < ",toString(Q3_pval)),
                                                                                      paste0("|logFC| >= ",toString(Q3_pval))),
                                                                             labels=c(paste0("|logFC| < ",toString(Q1_pval)," (",below_Q1_pval," DEGs)"),
                                                                                      paste0(toString(Q1_pval)," <= |logFC| < ",toString(Q2_pval)," (",below_Q1_Q2_pval," DEGs)"),
                                                                                      paste0(toString(Q2_pval)," <= |logFC| < ",toString(Q3_pval)," (",below_Q2_Q3_pval," DEGs)"),
                                                                                      paste0("|logFC| >= ",toString(Q3_pval)," (",above_Q3_pval," DEGs)")))))+
            geom_boxplot()+
            ggtitle("False Discovery Rate (FDR) when down-sampling samples at a range of values, by effect sizes (p-val)")+
            labs(y = "FDR", x = "Number of samples", fill="Effect size range") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9))+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "FPs_effects_pval.png"),FPs_effects_pval.plot,width=30,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "FPs_effects_pval.pdf"),FPs_effects_pval.plot,width=30,height=25,units="cm",bg="white")

        ## save all FDR/pval plots as 2 subplots
        # DEGs detected
        ggsave(file.path(savepath, "DEGs_detected_boxplot.png"),arrangeGrob(degs_detected_boxplot_fdr.plot,degs_detected_boxplot_pval.plot,ncol=2),width=50,height=20,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_detected_boxplot.pdf"),arrangeGrob(degs_detected_boxplot_fdr.plot,degs_detected_boxplot_pval.plot,ncol=2),width=50,height=20,units="cm",bg="white")
        # DEGs detected effect sizes
        ggsave(file.path(savepath, "DEGs_effects.png"),arrangeGrob(degs_effects_fdr.plot,degs_effects_pval.plot,ncol=2),width=70,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_effects.pdf"),arrangeGrob(degs_effects_fdr.plot,degs_effects_pval.plot,ncol=2),width=70,height=25,units="cm",bg="white")
        # upreg DEGs detected effect sizes
        ggsave(file.path(savepath, "DEGs_upreg_effects.png"),arrangeGrob(degs_effects_fdr_upreg.plot,degs_effects_pval_upreg.plot,ncol=2),width=70,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_upreg_effects.pdf"),arrangeGrob(degs_effects_fdr_upreg.plot,degs_effects_pval_upreg.plot,ncol=2),width=70,height=25,units="cm",bg="white")
        # downreg DEGs detected effect sizes
        ggsave(file.path(savepath, "DEGs_downreg_effects.png"),arrangeGrob(degs_effects_fdr_downreg.plot,degs_effects_pval_downreg.plot,ncol=2),width=70,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_downreg_effects.pdf"),arrangeGrob(degs_effects_fdr_downreg.plot,degs_effects_pval_downreg.plot,ncol=2),width=70,height=25,units="cm",bg="white")
        # FPs detected
        ggsave(file.path(savepath, "FPs_detected_boxplot.png"),arrangeGrob(FPs_detected_boxplot_fdr.plot,FPs_detected_boxplot_pval.plot,ncol=2),width=50,height=20,units="cm",bg="white")
        ggsave(file.path(savepath, "FPs_detected_boxplot.pdf"),arrangeGrob(FPs_detected_boxplot_fdr.plot,FPs_detected_boxplot_pval.plot,ncol=2),width=50,height=20,units="cm",bg="white")
        # FDR
        ggsave(file.path(savepath, "FDRs_boxplot.png"),arrangeGrob(FDRs_boxplot_fdr.plot,FDRs_boxplot_pval.plot,ncol=2),width=50,height=20,units="cm",bg="white")
        ggsave(file.path(savepath, "FDRs_boxplot.pdf"),arrangeGrob(FDRs_boxplot_fdr.plot,FDRs_boxplot_pval.plot,ncol=2),width=50,height=20,units="cm",bg="white")
        # FDR by effect size
        ggsave(file.path(savepath, "FPs_effects.png"),arrangeGrob(FPs_effects_fdr.plot,FPs_effects_pval.plot,ncol=2),width=60,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "FPs_effects.pdf"),arrangeGrob(FPs_effects_fdr.plot,FPs_effects_pval.plot,ncol=2),width=60,height=25,units="cm",bg="white")

    }else{

        #### plots using DEGs selected with an FDR cut-off
        # define path and folders        
        path <- file.path(inpath,"DE_downsampling_cells/")
        savepath <- file.path(output_path,"DE_downsampling_cells/")
        dir.create(savepath,showWarnings=FALSE,recursive=TRUE)
        downsampled_folders <- paste0(paste(range_downsampled, sep=" ", collapse=NULL),"cells_persample")
        ## plot % DEGs detected
        # load df with number of DEGs for each iteration/number of samples
        load(file.path(path, "DEGs_detected_fdr.RData"))
        # remove iteration column
        DEGs_detected_fdr <- subset(DEGs_detected_fdr,select=-c(Iteration))
        # convert to percentage
        DEGs_detected_fdr <- round(DEGs_detected_fdr/totalDEGs_fdr * 100, 2)
        # rename columns
        colnames(DEGs_detected_fdr) <- range_downsampled
        # reshape data
        DEGs_detected_new_fdr <- melt(DEGs_detected_fdr)
        ## create and save boxplots plot
        degs_detected_boxplot_cells_fdr.plot<-
            ggplot(data=DEGs_detected_new_fdr,
                    aes(x = factor(variable), y = value))+
            geom_boxplot(fill="#69b3a2",color="black")+
            ggtitle("Percentage of total DEGs detected when down-sampling cells at a range of values (FDR)")+
            labs(y = "% DEGs detected", x = "Mean number of cells per sample", fill="Cells") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9),legend.position="none")+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "DEGs_detected_boxplot_cells_fdr.png"),degs_detected_boxplot_cells_fdr.plot,width=30,height=20,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_detected_boxplot_cells_fdr.pdf"),degs_detected_boxplot_cells_fdr.plot,width=30,height=20,units="cm",bg="white")

        ## plot % of DEGs detected in each effect size range
        setwd(path)
        # loop through directories
        allDEGs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                DEGcounts <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], adj_pval<fdr)
                    # add numDEGs for each permutation
                    DEGcounts <- unlist(c(DEGcounts,
                                        length(subset(degs,name%in%DEGs_below_Q1_fdr)$name),
                                        length(subset(degs,name%in%DEGs_between_Q1_Q2_fdr)$name),
                                        length(subset(degs,name%in%DEGs_between_Q2_Q3_fdr)$name),
                                        length(subset(degs,name%in%DEGs_above_Q3_fdr)$name)))
                }
                allDEGs <- unlist(c(allDEGs, DEGcounts))
            }
        }
        # split by number of samples
        Iteration <- c()
        for(perm in 1:Nperms){
            Iteration <- c(Iteration, paste0(toString(perm),"_1"),
                                      paste0(toString(perm),"_2"),
                                      paste0(toString(perm),"_3"),
                                      paste0(toString(perm),"_4"))
        }
        # create dataframe
        allDEGs_effects <- data.frame(allDEGs[1:(4*Nperms)])
        colnames(allDEGs_effects) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            allDEGs_effects[[toString(range_downsampled[[value]])]] <- allDEGs[(4*(value-1)*Nperms+1):(4*value*Nperms)]
        }
        # convert to percentages
        for(j in 1:dim(allDEGs_effects)[[1]]){
            if(j%%4==1){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_fdr*100,2)
            }else if(j%%4==2){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_Q2_fdr*100,2)
            }else if(j%%4==3){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q2_Q3_fdr*100,2)
            }else{
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/above_Q3_fdr*100,2)
            }
        }
        # add and change iteration column above
        allDEGs_effects$`Iteration` <- Iteration
        allDEGs_effects$Iteration <- str_sub(allDEGs_effects$Iteration,-1,-1)
        # change name of iteration to be range of LogFCs
        allDEGs_effects$Iteration <- gsub("^1$",paste0("|logFC| < ",toString(Q1_fdr)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^2$",paste0(toString(Q1_fdr)," <= |logFC| < ",toString(Q2_fdr)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^3$",paste0(toString(Q2_fdr)," <= |logFC| < ",toString(Q3_fdr)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^4$",paste0("|logFC| >= ",toString(Q3_fdr)),allDEGs_effects$Iteration)
        # reshape
        allDEGs_effects_new <- melt(allDEGs_effects)
        # plot
        degs_cells_effects_fdr.plot <- 
            ggplot(data=allDEGs_effects_new,
                aes(x = factor(variable), y = value, fill = factor(Iteration,levels=c(paste0("|logFC| < ",toString(Q1_fdr)),
                                                                                      paste0(toString(Q1_fdr)," <= |logFC| < ",toString(Q2_fdr)),
                                                                                      paste0(toString(Q2_fdr)," <= |logFC| < ",toString(Q3_fdr)),
                                                                                      paste0("|logFC| >= ",toString(Q3_fdr))),
                                                                             labels=c(paste0("|logFC| < ",toString(Q1_fdr)," (",below_Q1_fdr," DEGs)"),
                                                                                      paste0(toString(Q1_fdr)," <= |logFC| < ",toString(Q2_fdr)," (",below_Q1_Q2_fdr," DEGs)"),
                                                                                      paste0(toString(Q2_fdr)," <= |logFC| < ",toString(Q3_fdr)," (",below_Q2_Q3_fdr," DEGs)"),
                                                                                      paste0("|logFC| >= ",toString(Q3_fdr)," (",above_Q3_fdr," DEGs)")))))+
            geom_boxplot()+
            ggtitle("Percentage of total DEGs detected when down-sampling cells at a range of values, by effect sizes (FDR)")+
            labs(y = "% DEGs detected", x = "Mean number of cells per sample", fill="Effect size range") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9))+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "DEGs_cells_effects_fdr.png"),degs_cells_effects_fdr.plot,width=30,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_cells_effects_fdr.pdf"),degs_cells_effects_fdr.plot,width=30,height=25,units="cm",bg="white")

        ## upreg plot
        setwd(path)
        # loop through directories
        allDEGs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                DEGcounts <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], adj_pval<fdr&logFC>0)
                    # add numDEGs for each permutation
                    DEGcounts <- unlist(c(DEGcounts,
                                        length(subset(degs,name%in%DEGs_below_Q1_fdr_upreg)$name),
                                        length(subset(degs,name%in%DEGs_between_Q1_Q2_fdr_upreg)$name),
                                        length(subset(degs,name%in%DEGs_between_Q2_Q3_fdr_upreg)$name),
                                        length(subset(degs,name%in%DEGs_above_Q3_fdr_upreg)$name)))
                }
                allDEGs <- unlist(c(allDEGs, DEGcounts))
            }
        }
        # split by number of samples
        Iteration <- c()
        for(perm in 1:Nperms){
            Iteration <- c(Iteration, paste0(toString(perm),"_1"),
                                    paste0(toString(perm),"_2"),
                                    paste0(toString(perm),"_3"),
                                    paste0(toString(perm),"_4"))
        }
        # create dataframe
        allDEGs_effects <- data.frame(allDEGs[1:(4*Nperms)])
        colnames(allDEGs_effects) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            allDEGs_effects[[toString(range_downsampled[[value]])]] <- allDEGs[(4*(value-1)*Nperms+1):(4*value*Nperms)]
        }
        # convert to percentages
        for(j in 1:dim(allDEGs_effects)[[1]]){
            if(j%%4==1){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_fdr_upreg*100,2)
            }else if(j%%4==2){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_Q2_fdr_upreg*100,2)
            }else if(j%%4==3){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q2_Q3_fdr_upreg*100,2)
            }else{
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/above_Q3_fdr_upreg*100,2)
            }
        }
        # add and change iteration column above
        allDEGs_effects$`Iteration` <- Iteration
        allDEGs_effects$Iteration <- str_sub(allDEGs_effects$Iteration,-1,-1)
        # change name of iteration to be range of LogFCs
        allDEGs_effects$Iteration <- gsub("^1$",paste0("0 < logFC < ",toString(Q1_fdr_upreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^2$",paste0(toString(Q1_fdr_upreg)," <= logFC < ",toString(Q2_fdr_upreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^3$",paste0(toString(Q2_fdr_upreg)," <= logFC < ",toString(Q3_fdr_upreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^4$",paste0("logFC >= ",toString(Q3_fdr_upreg)),allDEGs_effects$Iteration)
        # reshape
        allDEGs_effects_new <- melt(allDEGs_effects)
        # plot
        degs_cells_effects_fdr_upreg.plot <- 
            ggplot(data=allDEGs_effects_new,
                aes(x = factor(variable), y = value, fill = factor(Iteration,levels=c(paste0("0 < logFC < ",toString(Q1_fdr_upreg)),
                                                                                        paste0(toString(Q1_fdr_upreg)," <= logFC < ",toString(Q2_fdr_upreg)),
                                                                                        paste0(toString(Q2_fdr_upreg)," <= logFC < ",toString(Q3_fdr_upreg)),
                                                                                        paste0("logFC >= ",toString(Q3_fdr_upreg))),
                                                                                labels=c(paste0("0 < logFC < ",toString(Q1_fdr_upreg)," (",below_Q1_fdr_upreg," DEGs)"),
                                                                                        paste0(toString(Q1_fdr_upreg)," <= logFC < ",toString(Q2_fdr_upreg)," (",below_Q1_Q2_fdr_upreg," DEGs)"),
                                                                                        paste0(toString(Q2_fdr_upreg)," <= logFC < ",toString(Q3_fdr_upreg)," (",below_Q2_Q3_fdr_upreg," DEGs)"),
                                                                                        paste0("logFC >= ",toString(Q3_fdr_upreg)," (",above_Q3_fdr_upreg," DEGs)")))))+
            geom_boxplot()+
            ggtitle("Percentage of total (up-reg) DEGs detected when down-sampling cells at a range of values, by effect sizes (FDR)")+
            labs(y = "% DEGs detected", x = "Mean number of cells per sample", fill="Effect size range") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9))+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "DEGs_cells_upreg_effects_fdr.png"),degs_cells_effects_fdr_upreg.plot,width=35,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_cells_upreg_effects_fdr.pdf"),degs_cells_effects_fdr_upreg.plot,width=35,height=25,units="cm",bg="white")
        
        ## downreg plot
        setwd(path)
        # loop through directories
        allDEGs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                DEGcounts <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], adj_pval<fdr&logFC<0)
                    # add numDEGs for each permutation
                    DEGcounts <- unlist(c(DEGcounts,
                                        length(subset(degs,name%in%DEGs_below_Q1_fdr_downreg)$name),
                                        length(subset(degs,name%in%DEGs_between_Q1_Q2_fdr_downreg)$name),
                                        length(subset(degs,name%in%DEGs_between_Q2_Q3_fdr_downreg)$name),
                                        length(subset(degs,name%in%DEGs_above_Q3_fdr_downreg)$name)))
                }
                allDEGs <- unlist(c(allDEGs, DEGcounts))
            }
        }
        # split by number of samples
        Iteration <- c()
        for(perm in 1:Nperms){
            Iteration <- c(Iteration, paste0(toString(perm),"_1"),
                                    paste0(toString(perm),"_2"),
                                    paste0(toString(perm),"_3"),
                                    paste0(toString(perm),"_4"))
        }
        # create dataframe
        allDEGs_effects <- data.frame(allDEGs[1:(4*Nperms)])
        colnames(allDEGs_effects) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            allDEGs_effects[[toString(range_downsampled[[value]])]] <- allDEGs[(4*(value-1)*Nperms+1):(4*value*Nperms)]
        }
        # convert to percentages
        for(j in 1:dim(allDEGs_effects)[[1]]){
            if(j%%4==1){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_fdr_downreg*100,2)
            }else if(j%%4==2){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_Q2_fdr_downreg*100,2)
            }else if(j%%4==3){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q2_Q3_fdr_downreg*100,2)
            }else{
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/above_Q3_fdr_downreg*100,2)
            }
        }
        # add and change iteration column above
        allDEGs_effects$`Iteration` <- Iteration
        allDEGs_effects$Iteration <- str_sub(allDEGs_effects$Iteration,-1,-1)
        # change name of iteration to be range of LogFCs
        allDEGs_effects$Iteration <- gsub("^1$",paste0("logFC < ",toString(Q1_fdr_downreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^2$",paste0(toString(Q1_fdr_downreg)," <= logFC < ",toString(Q2_fdr_downreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^3$",paste0(toString(Q2_fdr_downreg)," <= logFC < ",toString(Q3_fdr_downreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^4$",paste0("0 > logFC >= ",toString(Q3_fdr_downreg)),allDEGs_effects$Iteration)
        # reshape
        allDEGs_effects_new <- melt(allDEGs_effects)
        # plot
        degs_cells_effects_fdr_downreg.plot <- 
            ggplot(data=allDEGs_effects_new,
                aes(x = factor(variable), y = value, fill = factor(Iteration,levels=c(paste0("logFC < ",toString(Q1_fdr_downreg)),
                                                                                        paste0(toString(Q1_fdr_downreg)," <= logFC < ",toString(Q2_fdr_downreg)),
                                                                                        paste0(toString(Q2_fdr_downreg)," <= logFC < ",toString(Q3_fdr_downreg)),
                                                                                        paste0("0 > logFC >= ",toString(Q3_fdr_downreg))),
                                                                                labels=c(paste0("logFC < ",toString(Q1_fdr_downreg)," (",below_Q1_fdr_downreg," DEGs)"),
                                                                                        paste0(toString(Q1_fdr_downreg)," <= logFC < ",toString(Q2_fdr_downreg)," (",below_Q1_Q2_fdr_downreg," DEGs)"),
                                                                                        paste0(toString(Q2_fdr_downreg)," <= logFC < ",toString(Q3_fdr_downreg)," (",below_Q2_Q3_fdr_downreg," DEGs)"),
                                                                                        paste0("0 > logFC >= ",toString(Q3_fdr_downreg)," (",above_Q3_fdr_downreg," DEGs)")))))+
            geom_boxplot()+
            ggtitle("Percentage of total (down-reg) DEGs detected when down-sampling cells at a range of values, by effect sizes (FDR)")+
            labs(y = "% DEGs detected", x = "Mean number of cells per sample", fill="Effect size range") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9))+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "DEGs_cells_downreg_effects_fdr.png"),degs_cells_effects_fdr_downreg.plot,width=35,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_cells_downreg_effects_fdr.pdf"),degs_cells_effects_fdr_downreg.plot,width=35,height=25,units="cm",bg="white")

        ## plot False Positive Rate
        setwd(path)
        # loop through directories
        allFPRs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                FPRs <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], adj_pval<fdr)$name
                    allgenes_tmp <- get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]]
                    non_degs <- subset(allgenes_tmp,!name%in%degs)
                    TNs <- length(subset(non_degs,name%in%nonDEGs_fdr$name)$name)
                    # get FPR
                    numFPs <- length(degs[degs %in% nonDEGs_fdr$name])
                    FPR <- numFPs/(numFPs+TNs)
                    # add numFPs for each permutation
                    FPRs <- unlist(c(FPRs, FPR))
                }
                allFPRs <- unlist(c(allFPRs, FPRs))
            }
        }
        # create dataframe
        all_FPRs <- data.frame(allFPRs[1:Nperms])
        colnames(all_FPRs) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            all_FPRs[[toString(range_downsampled[[value]])]] <- allFPRs[((value-1)*Nperms+1):(value*Nperms)]
        }
        # reshape data
        all_FPRs_new <- melt(all_FPRs)
        # create and save boxplots
        FPs_detected_cells_boxplot_fdr.plot<-
            ggplot(data=all_FPRs_new,
                    aes(x = factor(variable), y = value))+
            geom_boxplot(fill="#f37735",color="black")+
            ggtitle("False Positive Rate (FPR) when down-sampling cells at a range of values (FDR)")+
            labs(y = "FPR", x = "Mean number of cells per sample", fill="Cells") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9),legend.position="none")+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "FPs_detected_cells_boxplot_fdr.png"),FPs_detected_cells_boxplot_fdr.plot,width=25,height=20,units="cm",bg="white")
        ggsave(file.path(savepath, "FPs_detected_cells_boxplot_fdr.pdf"),FPs_detected_cells_boxplot_fdr.plot,width=25,height=20,units="cm",bg="white")

        ## plot False Discovery Rate
        setwd(path)
        # loop through directories
        allFDRs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                FDRs <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], adj_pval<fdr)$name
                    # get true positives
                    numTPs <- length(degs[degs %in% DEGs_fdr])
                    # get FDR
                    numFPs <- length(degs[degs %in% nonDEGs_fdr$name])
                    FDR <- numFPs/(numFPs+numTPs)
                    # add FDR for each permutation
                    FDRs <- unlist(c(FDRs, FDR))
                }
                allFDRs <- unlist(c(allFDRs, FDRs))
            }
        }
        # create dataframe
        all_FDRs <- data.frame(allFDRs[1:Nperms])
        colnames(all_FDRs) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            all_FDRs[[toString(range_downsampled[[value]])]] <- allFDRs[((value-1)*Nperms+1):(value*Nperms)]
        }
        # reshape data
        all_FDRs_new <- melt(all_FDRs)
        # create and save boxplots
        FDRs_cells_boxplot_fdr.plot<-
            ggplot(data=all_FDRs_new,
                    aes(x = factor(variable), y = value))+
            geom_boxplot(fill="#f37735",color="black")+
            ggtitle("False Discovery Rate (FDR) when down-sampling cells at a range of values (FDR)")+
            labs(y = "FDR", x = "Mean number of cells per sample", fill="Cells") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9),legend.position="none")+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "FDRs_cells_boxplot_fdr.png"),FDRs_cells_boxplot_fdr.plot,width=25,height=20,units="cm",bg="white")
        ggsave(file.path(savepath, "FDRs_cells_boxplot_fdr.pdf"),FDRs_cells_boxplot_fdr.plot,width=25,height=20,units="cm",bg="white")

        ## plot False Discovery Rate of DEGs in each effect size range
        setwd(path)
        # loop through directories
        allFDRs_effs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                FDRs <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], adj_pval<fdr)
                    # get numTPs for each range
                    TP_belowQ1 <- length(subset(degs,name%in%DEGs_below_Q1_fdr)$name)
                    TP_Q1_to_Q2 <- length(subset(degs,name%in%DEGs_between_Q1_Q2_fdr)$name)
                    TP_Q2_to_Q3 <- length(subset(degs,name%in%DEGs_between_Q2_Q3_fdr)$name)
                    TP_aboveQ3 <- length(subset(degs,name%in%DEGs_above_Q3_fdr)$name)
                    # get FPs for each range
                    FP_belowQ1 <- length(subset(degs,name%in%nonDEGs_belowQ1_fdr)$name)
                    FP_Q1_to_Q2 <- length(subset(degs,name%in%nonDEGs_Q1_to_Q2_fdr)$name)
                    FP_Q2_to_Q3 <- length(subset(degs,name%in%nonDEGs_Q2_to_Q3_fdr)$name)
                    FP_aboveQ3 <- length(subset(degs,name%in%nonDEGs_aboveQ3_fdr)$name)
                    # add numDEGs for each permutation
                    FDRs <- unlist(c(FDRs,
                                    FP_belowQ1/(FP_belowQ1+TP_belowQ1),
                                    FP_Q1_to_Q2/(FP_Q1_to_Q2+TP_Q1_to_Q2),
                                    FP_Q2_to_Q3/(FP_Q2_to_Q3+TP_Q2_to_Q3),
                                    FP_aboveQ3/(FP_aboveQ3+TP_aboveQ3)))
                }
                allFDRs_effs <- unlist(c(allFDRs_effs, FDRs))
            }
        }
        # split by number of samples
        Iteration <- c()
        for(perm in 1:Nperms){
            Iteration <- c(Iteration, paste0(toString(perm),"_1"),
                                      paste0(toString(perm),"_2"),
                                      paste0(toString(perm),"_3"),
                                      paste0(toString(perm),"_4"))
        }
        # create dataframe
        all_FDR_effs <- data.frame(allFDRs_effs[1:(4*Nperms)])
        colnames(all_FDR_effs) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            all_FDR_effs[[toString(range_downsampled[[value]])]] <- allFDRs_effs[(4*(value-1)*Nperms+1):(4*value*Nperms)]
        }
        ## create boxplots for various ranges of effect sizes
        # add and change iteration column above
        all_FDR_effs$`Iteration` <- Iteration
        all_FDR_effs$Iteration <- str_sub(all_FDR_effs$Iteration,-1,-1)
        # change name of iteration to be range of LogFCs
        all_FDR_effs$Iteration <- gsub("^1$",paste0("|logFC| < ",toString(Q1_fdr)),all_FDR_effs$Iteration)
        all_FDR_effs$Iteration <- gsub("^2$",paste0(toString(Q1_fdr)," <= |logFC| < ",toString(Q2_fdr)),all_FDR_effs$Iteration)
        all_FDR_effs$Iteration <- gsub("^3$",paste0(toString(Q2_fdr)," <= |logFC| < ",toString(Q3_fdr)),all_FDR_effs$Iteration)
        all_FDR_effs$Iteration <- gsub("^4$",paste0("|logFC| >= ",toString(Q3_fdr)),all_FDR_effs$Iteration)
        # reshape
        all_FDR_effs_new <- melt(all_FDR_effs)
        # plot
        FPs_cells_effects_fdr.plot <- 
            ggplot(data=all_FDR_effs_new,
                aes(x = factor(variable), y = value, fill = factor(Iteration,levels=c(paste0("|logFC| < ",toString(Q1_fdr)),
                                                                                      paste0(toString(Q1_fdr)," <= |logFC| < ",toString(Q2_fdr)),
                                                                                      paste0(toString(Q2_fdr)," <= |logFC| < ",toString(Q3_fdr)),
                                                                                      paste0("|logFC| >= ",toString(Q3_fdr))),
                                                                             labels=c(paste0("|logFC| < ",toString(Q1_fdr)," (",below_Q1_fdr," DEGs)"),
                                                                                      paste0(toString(Q1_fdr)," <= |logFC| < ",toString(Q2_fdr)," (",below_Q1_Q2_fdr," DEGs)"),
                                                                                      paste0(toString(Q2_fdr)," <= |logFC| < ",toString(Q3_fdr)," (",below_Q2_Q3_fdr," DEGs)"),
                                                                                      paste0("|logFC| >= ",toString(Q3_fdr)," (",above_Q3_fdr," DEGs)")))))+
            geom_boxplot()+
            ggtitle("False Discovery Rate (FDR) when down-sampling cells at a range of values, by effect sizes (FDR)")+
            labs(y = "FDR", x = "Mean number of cells per sample", fill="Effect size range") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9))+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "FPs_cells_effects_fdr.png"),FPs_cells_effects_fdr.plot,width=30,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "FPs_cells_effects_fdr.pdf"),FPs_cells_effects_fdr.plot,width=30,height=25,units="cm",bg="white")

        #### plots using DEGs selected with a nominal p-value cut-off
        ## plot % DEGs detected
        # load df with number of DEGs for each iteration/number of samples
        load(file.path(path, "DEGs_detected_pval.RData"))
        # remove iteration column
        DEGs_detected_pval <- subset(DEGs_detected_pval,select=-c(Iteration))
        # convert to percentage
        DEGs_detected_pval <- round(DEGs_detected_pval/totalDEGs_pval * 100, 2)
        # rename columns
        colnames(DEGs_detected_pval) <- range_downsampled
        # reshape data
        DEGs_detected_new_pval <- melt(DEGs_detected_pval)
        ## create and save boxplots plot
        degs_detected_boxplot_cells_pval.plot<-
            ggplot(data=DEGs_detected_new_pval,
                    aes(x = factor(variable), y = value))+
            geom_boxplot(fill="#69b3a2",color="black")+
            ggtitle("Percentage of total DEGs detected when down-sampling cells at a range of values (p-val)")+
            labs(y = "% DEGs detected", x = "Mean number of cells per sample", fill="Cells") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9),legend.position="none")+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "DEGs_detected_boxplot_cells_pval.png"),degs_detected_boxplot_cells_pval.plot,width=30,height=20,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_detected_boxplot_cells_pval.pdf"),degs_detected_boxplot_cells_pval.plot,width=30,height=20,units="cm",bg="white")

        ## plot % of DEGs detected in each effect size range
        setwd(path)
        # loop through directories
        allDEGs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                DEGcounts <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], PValue<nom_pval)
                    # add numDEGs for each permutation
                    DEGcounts <- unlist(c(DEGcounts,
                                        length(subset(degs,name%in%DEGs_below_Q1_pval)$name),
                                        length(subset(degs,name%in%DEGs_between_Q1_Q2_pval)$name),
                                        length(subset(degs,name%in%DEGs_between_Q2_Q3_pval)$name),
                                        length(subset(degs,name%in%DEGs_above_Q3_pval)$name)))
                }
                allDEGs <- unlist(c(allDEGs, DEGcounts))
            }
        }
        # split by number of samples
        Iteration <- c()
        for(perm in 1:Nperms){
            Iteration <- c(Iteration, paste0(toString(perm),"_1"),
                                      paste0(toString(perm),"_2"),
                                      paste0(toString(perm),"_3"),
                                      paste0(toString(perm),"_4"))
        }
        # create dataframe
        allDEGs_effects <- data.frame(allDEGs[1:(4*Nperms)])
        colnames(allDEGs_effects) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            allDEGs_effects[[toString(range_downsampled[[value]])]] <- allDEGs[(4*(value-1)*Nperms+1):(4*value*Nperms)]
        }
        # convert to percentages
        for(j in 1:dim(allDEGs_effects)[[1]]){
            if(j%%4==1){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_pval*100,2)
            }else if(j%%4==2){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_Q2_pval*100,2)
            }else if(j%%4==3){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q2_Q3_pval*100,2)
            }else{
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/above_Q3_pval*100,2)
            }
        }
        # add and change iteration column above
        allDEGs_effects$`Iteration` <- Iteration
        allDEGs_effects$Iteration <- str_sub(allDEGs_effects$Iteration,-1,-1)
        # change name of iteration to be range of LogFCs
        allDEGs_effects$Iteration <- gsub("^1$",paste0("|logFC| < ",toString(Q1_pval)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^2$",paste0(toString(Q1_pval)," <= |logFC| < ",toString(Q2_pval)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^3$",paste0(toString(Q2_pval)," <= |logFC| < ",toString(Q3_pval)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^4$",paste0("|logFC| >= ",toString(Q3_pval)),allDEGs_effects$Iteration)
        # reshape
        allDEGs_effects_new <- melt(allDEGs_effects)
        # plot
        degs_cells_effects_pval.plot <- 
            ggplot(data=allDEGs_effects_new,
                aes(x = factor(variable), y = value, fill = factor(Iteration,levels=c(paste0("|logFC| < ",toString(Q1_pval)),
                                                                                      paste0(toString(Q1_pval)," <= |logFC| < ",toString(Q2_pval)),
                                                                                      paste0(toString(Q2_pval)," <= |logFC| < ",toString(Q3_pval)),
                                                                                      paste0("|logFC| >= ",toString(Q3_pval))),
                                                                             labels=c(paste0("|logFC| < ",toString(Q1_pval)," (",below_Q1_pval," DEGs)"),
                                                                                      paste0(toString(Q1_pval)," <= |logFC| < ",toString(Q2_pval)," (",below_Q1_Q2_pval," DEGs)"),
                                                                                      paste0(toString(Q2_pval)," <= |logFC| < ",toString(Q3_pval)," (",below_Q2_Q3_pval," DEGs)"),
                                                                                      paste0("|logFC| >= ",toString(Q3_pval)," (",above_Q3_pval," DEGs)")))))+
            geom_boxplot()+
            ggtitle("Percentage of total DEGs detected when down-sampling cells at a range of values, by effect sizes (p-val)")+
            labs(y = "% DEGs detected", x = "Mean number of cells per sample", fill="Effect size range") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9))+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "DEGs_cells_effects_pval.png"),degs_cells_effects_pval.plot,width=30,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_cells_effects_pval.pdf"),degs_cells_effects_pval.plot,width=30,height=25,units="cm",bg="white")

        ## upreg plot
        setwd(path)
        # loop through directories
        allDEGs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                DEGcounts <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], PValue<nom_pval&logFC>0)
                    # add numDEGs for each permutation
                    DEGcounts <- unlist(c(DEGcounts,
                                        length(subset(degs,name%in%DEGs_below_Q1_pval_upreg)$name),
                                        length(subset(degs,name%in%DEGs_between_Q1_Q2_pval_upreg)$name),
                                        length(subset(degs,name%in%DEGs_between_Q2_Q3_pval_upreg)$name),
                                        length(subset(degs,name%in%DEGs_above_Q3_pval_upreg)$name)))
                }
                allDEGs <- unlist(c(allDEGs, DEGcounts))
            }
        }
        # split by number of samples
        Iteration <- c()
        for(perm in 1:Nperms){
            Iteration <- c(Iteration, paste0(toString(perm),"_1"),
                                    paste0(toString(perm),"_2"),
                                    paste0(toString(perm),"_3"),
                                    paste0(toString(perm),"_4"))
        }
        # create dataframe
        allDEGs_effects <- data.frame(allDEGs[1:(4*Nperms)])
        colnames(allDEGs_effects) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            allDEGs_effects[[toString(range_downsampled[[value]])]] <- allDEGs[(4*(value-1)*Nperms+1):(4*value*Nperms)]
        }
        # convert to percentages
        for(j in 1:dim(allDEGs_effects)[[1]]){
            if(j%%4==1){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_pval_upreg*100,2)
            }else if(j%%4==2){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_Q2_pval_upreg*100,2)
            }else if(j%%4==3){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q2_Q3_pval_upreg*100,2)
            }else{
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/above_Q3_pval_upreg*100,2)
            }
        }
        # add and change iteration column above
        allDEGs_effects$`Iteration` <- Iteration
        allDEGs_effects$Iteration <- str_sub(allDEGs_effects$Iteration,-1,-1)
        # change name of iteration to be range of LogFCs
        allDEGs_effects$Iteration <- gsub("^1$",paste0("0 < logFC < ",toString(Q1_pval_upreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^2$",paste0(toString(Q1_pval_upreg)," <= logFC < ",toString(Q2_pval_upreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^3$",paste0(toString(Q2_pval_upreg)," <= logFC < ",toString(Q3_pval_upreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^4$",paste0("logFC >= ",toString(Q3_pval_upreg)),allDEGs_effects$Iteration)
        # reshape
        allDEGs_effects_new <- melt(allDEGs_effects)
        # plot
        degs_cells_effects_pval_upreg.plot <- 
            ggplot(data=allDEGs_effects_new,
                aes(x = factor(variable), y = value, fill = factor(Iteration,levels=c(paste0("0 < logFC < ",toString(Q1_pval_upreg)),
                                                                                        paste0(toString(Q1_pval_upreg)," <= logFC < ",toString(Q2_pval_upreg)),
                                                                                        paste0(toString(Q2_pval_upreg)," <= logFC < ",toString(Q3_pval_upreg)),
                                                                                        paste0("logFC >= ",toString(Q3_pval_upreg))),
                                                                                labels=c(paste0("0 < logFC < ",toString(Q1_pval_upreg)," (",below_Q1_pval_upreg," DEGs)"),
                                                                                        paste0(toString(Q1_pval_upreg)," <= logFC < ",toString(Q2_pval_upreg)," (",below_Q1_Q2_pval_upreg," DEGs)"),
                                                                                        paste0(toString(Q2_pval_upreg)," <= logFC < ",toString(Q3_pval_upreg)," (",below_Q2_Q3_pval_upreg," DEGs)"),
                                                                                        paste0("logFC >= ",toString(Q3_pval_upreg)," (",above_Q3_pval_upreg," DEGs)")))))+
            geom_boxplot()+
            ggtitle("Percentage of total (up-reg) DEGs detected when down-sampling cells at a range of values, by effect sizes (p-val)")+
            labs(y = "% DEGs detected", x = "Mean number of cells per sample", fill="Effect size range") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9))+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "DEGs_cells_upreg_effects_pval.png"),degs_cells_effects_pval_upreg.plot,width=35,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_cells_upreg_effects_pval.pdf"),degs_cells_effects_pval_upreg.plot,width=35,height=25,units="cm",bg="white")
        
        ## downreg plot
        setwd(path)
        # loop through directories
        allDEGs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                DEGcounts <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], PValue<nom_pval&logFC<0)
                    # add numDEGs for each permutation
                    DEGcounts <- unlist(c(DEGcounts,
                                        length(subset(degs,name%in%DEGs_below_Q1_pval_downreg)$name),
                                        length(subset(degs,name%in%DEGs_between_Q1_Q2_pval_downreg)$name),
                                        length(subset(degs,name%in%DEGs_between_Q2_Q3_pval_downreg)$name),
                                        length(subset(degs,name%in%DEGs_above_Q3_pval_downreg)$name)))
                }
                allDEGs <- unlist(c(allDEGs, DEGcounts))
            }
        }
        # split by number of samples
        Iteration <- c()
        for(perm in 1:Nperms){
            Iteration <- c(Iteration, paste0(toString(perm),"_1"),
                                    paste0(toString(perm),"_2"),
                                    paste0(toString(perm),"_3"),
                                    paste0(toString(perm),"_4"))
        }
        # create dataframe
        allDEGs_effects <- data.frame(allDEGs[1:(4*Nperms)])
        colnames(allDEGs_effects) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            allDEGs_effects[[toString(range_downsampled[[value]])]] <- allDEGs[(4*(value-1)*Nperms+1):(4*value*Nperms)]
        }
        # convert to percentages
        for(j in 1:dim(allDEGs_effects)[[1]]){
            if(j%%4==1){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_pval_downreg*100,2)
            }else if(j%%4==2){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q1_Q2_pval_downreg*100,2)
            }else if(j%%4==3){
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/below_Q2_Q3_pval_downreg*100,2)
            }else{
                allDEGs_effects[j,] <- round(allDEGs_effects[j,]/above_Q3_pval_downreg*100,2)
            }
        }
        # add and change iteration column above
        allDEGs_effects$`Iteration` <- Iteration
        allDEGs_effects$Iteration <- str_sub(allDEGs_effects$Iteration,-1,-1)
        # change name of iteration to be range of LogFCs
        allDEGs_effects$Iteration <- gsub("^1$",paste0("logFC < ",toString(Q1_pval_downreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^2$",paste0(toString(Q1_pval_downreg)," <= logFC < ",toString(Q2_pval_downreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^3$",paste0(toString(Q2_pval_downreg)," <= logFC < ",toString(Q3_pval_downreg)),allDEGs_effects$Iteration)
        allDEGs_effects$Iteration <- gsub("^4$",paste0("0 > logFC >= ",toString(Q3_pval_downreg)),allDEGs_effects$Iteration)
        # reshape
        allDEGs_effects_new <- melt(allDEGs_effects)
        # plot
        degs_cells_effects_pval_downreg.plot <- 
            ggplot(data=allDEGs_effects_new,
                aes(x = factor(variable), y = value, fill = factor(Iteration,levels=c(paste0("logFC < ",toString(Q1_pval_downreg)),
                                                                                        paste0(toString(Q1_pval_downreg)," <= logFC < ",toString(Q2_pval_downreg)),
                                                                                        paste0(toString(Q2_pval_downreg)," <= logFC < ",toString(Q3_pval_downreg)),
                                                                                        paste0("0 > logFC >= ",toString(Q3_pval_downreg))),
                                                                                labels=c(paste0("logFC < ",toString(Q1_pval_downreg)," (",below_Q1_pval_downreg," DEGs)"),
                                                                                        paste0(toString(Q1_pval_downreg)," <= logFC < ",toString(Q2_pval_downreg)," (",below_Q1_Q2_pval_downreg," DEGs)"),
                                                                                        paste0(toString(Q2_pval_downreg)," <= logFC < ",toString(Q3_pval_downreg)," (",below_Q2_Q3_pval_downreg," DEGs)"),
                                                                                        paste0("0 > logFC >= ",toString(Q3_pval_downreg)," (",above_Q3_pval_downreg," DEGs)")))))+
            geom_boxplot()+
            ggtitle("Percentage of total (down-reg) DEGs detected when down-sampling cells at a range of values, by effect sizes (p-val)")+
            labs(y = "% DEGs detected", x = "Mean number of cells per sample", fill="Effect size range") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9))+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "DEGs_cells_downreg_effects_pval.png"),degs_cells_effects_pval_downreg.plot,width=35,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_cells_downreg_effects_pval.pdf"),degs_cells_effects_pval_downreg.plot,width=35,height=25,units="cm",bg="white")

        ## plot False Positive Rate
        setwd(path)
        # loop through directories
        allFPRs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                FPRs <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], PValue<nom_pval)$name
                    allgenes_tmp <- get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]]
                    non_degs <- subset(allgenes_tmp,!name%in%degs)
                    TNs <- length(subset(non_degs,name%in%nonDEGs_pval$name)$name)
                    # get FPR
                    numFPs <- length(degs[degs %in% nonDEGs_pval$name])
                    FPR <- numFPs/(numFPs+TNs)
                    # add numFPs for each permutation
                    FPRs <- unlist(c(FPRs, FPR))
                }
                allFPRs <- unlist(c(allFPRs, FPRs))
            }
        }
        # create dataframe
        all_FPRs <- data.frame(allFPRs[1:Nperms])
        colnames(all_FPRs) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            all_FPRs[[toString(range_downsampled[[value]])]] <- allFPRs[((value-1)*Nperms+1):(value*Nperms)]
        }
        # reshape data
        all_FPRs_new <- melt(all_FPRs)
        # create and save boxplots
        FPs_detected_cells_boxplot_pval.plot<-
            ggplot(data=all_FPRs_new,
                    aes(x = factor(variable), y = value))+
            geom_boxplot(fill="#f37735",color="black")+
            ggtitle("False Positive Rate (FPR) when down-sampling cells at a range of values (p-val)")+
            labs(y = "FPR", x = "Mean number of cells per sample", fill="Cells") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9),legend.position="none")+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "FPs_detected_cells_boxplot_pval.png"),FPs_detected_cells_boxplot_pval.plot,width=25,height=20,units="cm",bg="white")
        ggsave(file.path(savepath, "FPs_detected_cells_boxplot_pval.pdf"),FPs_detected_cells_boxplot_pval.plot,width=25,height=20,units="cm",bg="white")

        ## plot False Discovery Rate
        setwd(path)
        # loop through directories
        allFDRs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                FDRs <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], PValue<nom_pval)$name
                    # get true positives
                    numTPs <- length(degs[degs %in% DEGs_pval])
                    # get FDR
                    numFPs <- length(degs[degs %in% nonDEGs_pval$name])
                    FDR <- numFPs/(numFPs+numTPs)
                    # add FDR for each permutation
                    FDRs <- unlist(c(FDRs, FDR))
                }
                allFDRs <- unlist(c(allFDRs, FDRs))
            }
        }
        # create dataframe
        all_FDRs <- data.frame(allFDRs[1:Nperms])
        colnames(all_FDRs) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            all_FDRs[[toString(range_downsampled[[value]])]] <- allFDRs[((value-1)*Nperms+1):(value*Nperms)]
        }
        # reshape data
        all_FDRs_new <- melt(all_FDRs)
        # create and save boxplots
        FDRs_cells_boxplot_pval.plot<-
            ggplot(data=all_FDRs_new,
                    aes(x = factor(variable), y = value))+
            geom_boxplot(fill="#f37735",color="black")+
            ggtitle("False Discovery Rate (FDR) when down-sampling cells at a range of values (p-val)")+
            labs(y = "FDR", x = "Mean number of cells per sample", fill="Cells") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9),legend.position="none")+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "FDRs_cells_boxplot_pval.png"),FDRs_cells_boxplot_pval.plot,width=25,height=20,units="cm",bg="white")
        ggsave(file.path(savepath, "FDRs_cells_boxplot_pval.pdf"),FDRs_cells_boxplot_pval.plot,width=25,height=20,units="cm",bg="white")

        ## plot False Discovery Rate of DEGs in each effect size range
        setwd(path)
        # loop through directories
        allFDRs_effs <- list()
        for(folder in mixedsort(list.files())){
            # remove other files above
            if(folder%in%downsampled_folders){
                newpath <- file.path(path,folder)
                # enter directory
                setwd(newpath)
                FDRs <- list()
                for(subfolder in mixedsort(list.files())){
                    print(paste0(folder,",",subfolder))
                    subpath <- file.path(newpath,subfolder)
                    # enter subdirectory
                    setwd(subpath)
                    # get samples
                    samples <- strsplit(subfolder,"_")[[1]][[1]]
                    # read DGE analysis output
                    load(paste0("DEout",subfolder,".RData"))
                    degs <- subset(get(paste0("DEout_",samples))$celltype_all_genes[[celltype_name]], PValue<nom_pval)
                    # get numTPs for each range
                    TP_belowQ1 <- length(subset(degs,name%in%DEGs_below_Q1_pval)$name)
                    TP_Q1_to_Q2 <- length(subset(degs,name%in%DEGs_between_Q1_Q2_pval)$name)
                    TP_Q2_to_Q3 <- length(subset(degs,name%in%DEGs_between_Q2_Q3_pval)$name)
                    TP_aboveQ3 <- length(subset(degs,name%in%DEGs_above_Q3_pval)$name)
                    # get FPs for each range
                    FP_belowQ1 <- length(subset(degs,name%in%nonDEGs_belowQ1_pval)$name)
                    FP_Q1_to_Q2 <- length(subset(degs,name%in%nonDEGs_Q1_to_Q2_pval)$name)
                    FP_Q2_to_Q3 <- length(subset(degs,name%in%nonDEGs_Q2_to_Q3_pval)$name)
                    FP_aboveQ3 <- length(subset(degs,name%in%nonDEGs_aboveQ3_pval)$name)
                    # add numDEGs for each permutation
                    FDRs <- unlist(c(FDRs,
                                    FP_belowQ1/(FP_belowQ1+TP_belowQ1),
                                    FP_Q1_to_Q2/(FP_Q1_to_Q2+TP_Q1_to_Q2),
                                    FP_Q2_to_Q3/(FP_Q2_to_Q3+TP_Q2_to_Q3),
                                    FP_aboveQ3/(FP_aboveQ3+TP_aboveQ3)))
                }
                allFDRs_effs <- unlist(c(allFDRs_effs, FDRs))
            }
        }
        # split by number of samples
        Iteration <- c()
        for(perm in 1:Nperms){
            Iteration <- c(Iteration, paste0(toString(perm),"_1"),
                                      paste0(toString(perm),"_2"),
                                      paste0(toString(perm),"_3"),
                                      paste0(toString(perm),"_4"))
        }
        # create dataframe
        all_FDR_effs <- data.frame(allFDRs_effs[1:(4*Nperms)])
        colnames(all_FDR_effs) <- c(toString(range_downsampled[[1]]))
        # loop through all values downsampled at
        for(value in 2:length(range_downsampled)){
            all_FDR_effs[[toString(range_downsampled[[value]])]] <- allFDRs_effs[(4*(value-1)*Nperms+1):(4*value*Nperms)]
        }
        ## create boxplots for various ranges of effect sizes
        # add and change iteration column above
        all_FDR_effs$`Iteration` <- Iteration
        all_FDR_effs$Iteration <- str_sub(all_FDR_effs$Iteration,-1,-1)
        # change name of iteration to be range of LogFCs
        all_FDR_effs$Iteration <- gsub("^1$",paste0("|logFC| < ",toString(Q1_pval)),all_FDR_effs$Iteration)
        all_FDR_effs$Iteration <- gsub("^2$",paste0(toString(Q1_pval)," <= |logFC| < ",toString(Q2_pval)),all_FDR_effs$Iteration)
        all_FDR_effs$Iteration <- gsub("^3$",paste0(toString(Q2_pval)," <= |logFC| < ",toString(Q3_pval)),all_FDR_effs$Iteration)
        all_FDR_effs$Iteration <- gsub("^4$",paste0("|logFC| >= ",toString(Q3_pval)),all_FDR_effs$Iteration)
        # reshape
        all_FDR_effs_new <- melt(all_FDR_effs)
        # plot
        FPs_cells_effects_pval.plot <- 
            ggplot(data=all_FDR_effs_new,
                aes(x = factor(variable), y = value, fill = factor(Iteration,levels=c(paste0("|logFC| < ",toString(Q1_pval)),
                                                                                      paste0(toString(Q1_pval)," <= |logFC| < ",toString(Q2_pval)),
                                                                                      paste0(toString(Q2_pval)," <= |logFC| < ",toString(Q3_pval)),
                                                                                      paste0("|logFC| >= ",toString(Q3_pval))),
                                                                             labels=c(paste0("|logFC| < ",toString(Q1_pval)," (",below_Q1_pval," DEGs)"),
                                                                                      paste0(toString(Q1_pval)," <= |logFC| < ",toString(Q2_pval)," (",below_Q1_Q2_pval," DEGs)"),
                                                                                      paste0(toString(Q2_pval)," <= |logFC| < ",toString(Q3_pval)," (",below_Q2_Q3_pval," DEGs)"),
                                                                                      paste0("|logFC| >= ",toString(Q3_pval)," (",above_Q3_pval," DEGs)")))))+
            geom_boxplot()+
            ggtitle("False Discovery Rate (FDR) when down-sampling cells at a range of values, by effect sizes (p-val)")+
            labs(y = "FDR", x = "Mean number of cells per sample", fill="Effect size range") +
            theme_cowplot()+
            theme(axis.text = element_text(size=9))+
            scale_fill_viridis(discrete = T)
        ggsave(file.path(savepath, "FPs_cells_effects_pval.png"),FPs_cells_effects_pval.plot,width=30,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "FPs_cells_effects_pval.pdf"),FPs_cells_effects_pval.plot,width=30,height=25,units="cm",bg="white")

        ## save all FDR/pval plots as 2 subplots
        # DEGs detected
        ggsave(file.path(savepath, "DEGs_detected_boxplot_cells.png"),arrangeGrob(degs_detected_boxplot_cells_fdr.plot,degs_detected_boxplot_cells_pval.plot,ncol=2),width=60,height=20,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_detected_boxplot_cells.pdf"),arrangeGrob(degs_detected_boxplot_cells_fdr.plot,degs_detected_boxplot_cells_pval.plot,ncol=2),width=60,height=20,units="cm",bg="white")
        # DEGs detected effect sizes
        ggsave(file.path(savepath, "DEGs_cells_effects.png"),arrangeGrob(degs_cells_effects_fdr.plot,degs_cells_effects_pval.plot,ncol=2),width=60,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_cells_effects.pdf"),arrangeGrob(degs_cells_effects_fdr.plot,degs_cells_effects_pval.plot,ncol=2),width=60,height=25,units="cm",bg="white")
        # upreg DEGs detected effect sizes
        ggsave(file.path(savepath, "DEGs_cells_upreg_effects.png"),arrangeGrob(degs_cells_effects_fdr_upreg.plot,degs_cells_effects_pval_upreg.plot,ncol=2),width=70,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_cells_upreg_effects.pdf"),arrangeGrob(degs_cells_effects_fdr_upreg.plot,degs_cells_effects_pval_upreg.plot,ncol=2),width=70,height=25,units="cm",bg="white")
        # downreg DEGs detected effect sizes
        ggsave(file.path(savepath, "DEGs_cells_downreg_effects.png"),arrangeGrob(degs_cells_effects_fdr_downreg.plot,degs_cells_effects_pval_downreg.plot,ncol=2),width=70,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "DEGs_cells_downreg_effects.pdf"),arrangeGrob(degs_cells_effects_fdr_downreg.plot,degs_cells_effects_pval_downreg.plot,ncol=2),width=70,height=25,units="cm",bg="white")
        # FPs detected
        ggsave(file.path(savepath, "FPs_detected_cells_boxplot.png"),arrangeGrob(FPs_detected_cells_boxplot_fdr.plot,FPs_detected_cells_boxplot_pval.plot,ncol=2),width=50,height=20,units="cm",bg="white")
        ggsave(file.path(savepath, "FPs_detected_cells_boxplot.pdf"),arrangeGrob(FPs_detected_cells_boxplot_fdr.plot,FPs_detected_cells_boxplot_pval.plot,ncol=2),width=50,height=20,units="cm",bg="white")
        # FDR
        ggsave(file.path(savepath, "FDRs_cells_boxplot.png"),arrangeGrob(FDRs_cells_boxplot_fdr.plot,FDRs_cells_boxplot_pval.plot,ncol=2),width=50,height=20,units="cm",bg="white")
        ggsave(file.path(savepath, "FDRs_cells_boxplot.pdf"),arrangeGrob(FDRs_cells_boxplot_fdr.plot,FDRs_cells_boxplot_pval.plot,ncol=2),width=50,height=20,units="cm",bg="white")
        # FDR by effect size
        ggsave(file.path(savepath, "FPs_cells_effects.png"),arrangeGrob(FPs_cells_effects_fdr.plot,FPs_cells_effects_pval.plot,ncol=2),width=60,height=25,units="cm",bg="white")
        ggsave(file.path(savepath, "FPs_cells_effects.pdf"),arrangeGrob(FPs_cells_effects_fdr.plot,FPs_cells_effects_pval.plot,ncol=2),width=60,height=25,units="cm",bg="white")

    }

}