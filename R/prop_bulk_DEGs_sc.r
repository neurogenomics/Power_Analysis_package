# Define global variables
utils::globalVariables(c("DEGs","numSamples","pctDEGs"))

#' Obtain overall percentage overlap between DEGs from bulk data (DE across all tissues) and various scRNA-seq datasets, across all common cell types

#' @importFrom stringr str_split
#' @importFrom utils tail
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 ggplot labs theme aes element_text geom_boxplot position_dodge2 scale_y_continuous guides guide_legend scale_fill_manual scale_alpha
#' @importFrom cowplot theme_cowplot

#' @param bulkDE DGE analysis output for a bulk RNA-seq dataset (e.g., `LFSR.tsv`): rows (rownames) should be the genes, columns should be tissues, and entries should be significance levels
#' @param bulk_cutoff Numeric. Proportion (0–1) of bulk tissues in which a gene must be differentially expressed to be considered (e.g., 0.9 selects DEGs found in ≥90% of tissues).
#' @param pvalue Numeric. P-value threshold for defining DEGs in the bulk dataset.
#' @param sampled Specifies the unit of down-sampling. Can be either `"individuals"` or `"cells"`, depending on whether the analysis downsamples across samples or cells.
#' @param fontsize_axislabels font size for axis labels in plot
#' @param fontsize_axisticks font size for axis tick labels in plot
#' @param fontsize_title font size for plot title
#' @param fontsize_legendlabels font size for legend labels in plot
#' @param fontsize_legendtitle font size for legend title in plot
#' @param plot_title plot title
#' @param output_path A directory path where down-sampled outputs and diagnostic plots will be saved.

#' Saves plot showing percentage DEGs from bulk data found in each scRNA-seq dataset, in the appropriate directory

prop_bulk_DEGs_sc <- function(bulkDE,
                              bulk_cutoff=0.9,
                              pvalue=0.05,
                              sampled="individuals",
                              fontsize_axislabels=12,
                              fontsize_axisticks=9,
                              fontsize_title=14,
                              fontsize_legendlabels=9,
                              fontsize_legendtitle=9,
                              plot_title="placeholder",
                              output_path=getwd()){

    sampled <- match.arg(sampled, choices = c("individuals", "cells"))

    # validate function input params
    validate_input_parameters_bulk(bulkDE=bulkDE, output_path=output_path,
                                   bulk_cutoff=bulk_cutoff,pvalue=pvalue, fontsize_axislabels=fontsize_axislabels,
                                   fontsize_axisticks=fontsize_axisticks, fontsize_title=fontsize_title, fontsize_legendlabels=fontsize_legendlabels,
                                   fontsize_legendtitle=fontsize_legendtitle, plot_title=plot_title)

    # default placeholder
    if(plot_title=="placeholder"){
        plot_title="Percentage DEGs from Bulk data detected when down-sampling scRNA-seq datasets (across all cell types)"
    }

    ## get DEGs significant in bulk data across at least bulk_cutoff of the tissues
    # deal with NA rows
    bulkDE[is.na(bulkDE)]<-100
    # total number of tissues
    numTissues <- length(bulkDE)
    # get DEGs
    bulk_DEGs <- c()
    for(rownum in 1:nrow(bulkDE)){
        if(sum(1*(bulkDE[rownum,]<pvalue))>=round(bulk_cutoff*numTissues)){
            bulk_DEGs <- c(bulk_DEGs,rownames(bulkDE[rownum,]))
        }
    }
    # get DEG names
    bulk_DEGs <- sapply(strsplit(bulk_DEGs,".",1), `[`, 1)

    ## get proportion of DEGs
    # lists to hold each var
    Dataset <- c()
    NumSamples <- c()
    Perm <- c()
    PctDEGs <- c()
    # reference the appropriate subdirectory based on the down-sampling unit
    folder_tag <- if (sampled == "individuals") "DE_downsampling" else "DE_downsampling_cells"
    for(dataset in list.dirs(output_path,recursive=F,full.names=F)){
        print(paste0("Downsampling ",dataset))
        # go inside dataset directory overall folder, down-sampling folder
        data_dir <- file.path(output_path, dataset, "Overall", folder_tag)
        downsampled_folders <- list.dirs(data_dir, recursive=FALSE, full.names=FALSE)
        # go inside numSamples, if exists
        for(sample_samples in downsampled_folders){
            sample <- as.numeric(gsub("[^0-9]", "", sample_samples))
            print(paste0("Processing ", sample_samples, "..."))
            # perms
            data_sample_dir <- file.path(data_dir,sample_samples)
            # check if this sample point exists
            if(dir.exists(data_sample_dir)){
                # loop through perms
                for(sample_perm in list.dirs(data_sample_dir,recursive=F,full.names=F)){
                    data_sample_perm_dir <- file.path(data_sample_dir,sample_perm)
                    # load DGE analysis output for this permutation
                    load(file.path(data_sample_perm_dir,paste0("DEout",sample_perm,".RData")))
                    # get % DEGs from Bulk Data DE here, add to vars, name dataset appropriately
                    Dataset <- c(Dataset, dataset)
                    NumSamples <- c(NumSamples, sample)
                    Perm <- c(Perm, str_split(sample_perm,"_")[[1]][2])
                    PctDEGs <- c(PctDEGs, sum((bulk_DEGs%in%DEGs)*1)/length(bulk_DEGs))
                }
            }
        }
    }
    # put together df
    DEGs_df <- data.frame(Dataset)
    DEGs_df$numSamples <- NumSamples
    DEGs_df$perm <- Perm
    DEGs_df$pctDEGs <- PctDEGs

    # function to select colours appropriately for boxplots
    generate_color_palette <- function(N, palette = "Set1") {
        max_colors <- 9  # the maximum number of colors available for the chosen palette
        if (N <= max_colors) {
            return(brewer.pal(N, palette))
        }else{
            # set seed
            set.seed(123)
            # create a continuous range of colors from the chosen palette
            continuous_palette <- colorRampPalette(brewer.pal(max_colors, palette))
            # generate N colors from the continuous range
            return(continuous_palette(N))
        }
    }

    # specify x-axis label according to the down-sampling style
    x_label <- if (sampled == "individuals") {
        "Number of samples"
    } else {
        "Mean number of cells per sample"
    }

    # create and return plot
    propDEGs_overall.plot <- ggplot(DEGs_df,aes(x=factor(numSamples),y=100*pctDEGs))+
                                    geom_boxplot(outlier.shape=NA,aes(fill=factor(Dataset)),position = position_dodge2(width = 1, preserve = "single"))+
                                    scale_y_continuous(labels = function(x) round(x))+
                                    guides(fill = guide_legend(label.format = function(x) round(as.numeric(x))))+
                                    theme_cowplot()+
                                    scale_fill_manual(values=generate_color_palette(length(list.dirs(output_path,recursive=F,full.names=F)),palette="Set1"))+
                                    labs(y="% DEGs", x=x_label, fill="Dataset",title="Percentage DEGs from Bulk data detected when down-sampling scRNA-seq datasets (across all cell types)")+
                                    scale_alpha(guide = 'none')+
                                    theme(axis.title = element_text(size = fontsize_axislabels),
                                          axis.text.x = element_text(size = fontsize_axisticks),
                                          axis.text.y = element_text(size = fontsize_axisticks),
                                          plot.title = element_text(size = fontsize_title),
                                          legend.text = element_text(size = fontsize_legendlabels),
                                          legend.title = element_text(size = fontsize_legendtitle))

    # save plot
    ggsave(file.path(output_path,"prop_bulk_DEGs_sc.png"), propDEGs_overall.plot, width=40, height=15, units="cm", bg="white")
    ggsave(file.path(output_path,"prop_bulk_DEGs_sc.pdf"), propDEGs_overall.plot, width=40, height=15, units="cm", bg="white")

}
