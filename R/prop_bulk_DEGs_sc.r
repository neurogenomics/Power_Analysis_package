# Define global variables
utils::globalVariables(c("DEGs","numSamples","pctDEGs"))

#' Obtain overall percentage overlap between DEGs from bulk data (DE across all tissues) and various scRNA-seq datasets, across all cell types

#' @importFrom stringr str_split
#' @importFrom utils tail
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 ggplot labs theme aes element_text geom_boxplot position_dodge2 scale_y_continuous guides guide_legend scale_fill_manual scale_alpha
#' @importFrom cowplot theme_cowplot

#' @param bulkDE DGE analysis output for a bulk RNA-seq dataset: rows (rownames) should be the genes, columns should be tissues, and entries should be significance levels
#' @param inpath path storing the down-sampled DGE analysis outputs for each dataset
#' @param range_downsampled vector or list containing values which the data will be downsampled at, in ascending order
#' @param bulk_cutoff percentage (proportion), specified so that we select DEGs common across >=bulk_cutoff of the tissues in the Bulk dataset
#' @param pvalue the cut-off p-value used to select DEGs (for bulk data)
#' @param fontsize_axislabels font size for axis labels in plot
#' @param fontsize_axisticks font size for axis tick labels in plot
#' @param fontsize_title font size for plot title
#' @param fontsize_legendlabels font size for legend labels in plot
#' @param fontsize_legendtitle font size for legend title in plot
#' @param plot_title plot title

#' @return plot showing percentage DEGs from bulk data found in each scRNA-seq dataset

prop_bulk_DEGs_sc <- function(bulkDE,
                              inpath,
                              range_downsampled,
                              bulk_cutoff=0.9,
                              pvalue=0.05,
                              fontsize_axislabels=12,
                              fontsize_axisticks=9,
                              fontsize_title=14,
                              fontsize_legendlabels=9,
                              fontsize_legendtitle=9,
                              plot_title="placeholder"){

    # validate function input params
    validate_input_parameters_bulk(bulkDE=bulkDE, inpath=inpath, range_downsampled=range_downsampled,
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
    for(dataset in list.dirs(inpath,recursive=F,full.names=F)){
        print(paste0("Downsampling ",dataset))
        # go inside dataset directory overall folder, downsampling folder
        data_dir <- paste0(inpath, "/", dataset, "/Overall/DE_downsampling/")
        setwd(data_dir)
        # go inside numSamples, if exists
        for(sample in range_downsampled){
            print(paste0(sample," samples"))
            # add "samples"
            sample_samples <- paste0(sample,"samples")
            data_sample_dir <- paste0(data_dir,sample_samples,"/")
            # check if this sample point exists
            if(dir.exists(data_sample_dir)){
                # go into directory
                setwd(data_sample_dir)
                # loop through perms
                for(perm in list.dirs(data_sample_dir,recursive=F,full.names=F)){
                    data_sample_perm_dir <- paste0(data_sample_dir,perm)
                    # go into each perm
                    setwd(data_sample_perm_dir)
                    # load DGE analysis output for this permutation
                    load(paste0("DEout",perm,".RData"))
                    # get % DEGs from Bulk Data DE here, add to vars, name dataset appropriately
                    if(length(str_split(dataset,"_")[[1]])==1){
                        dataset_name <- sub("^\\w", toupper(substring(dataset, 1, 1)), dataset)
                    }else{
                        split_name <- str_split(sub("^\\w", toupper(substring(dataset, 1, 1)), dataset),"_")[[1]]
                        tmp_splitname <- split_name[-length(split_name)]
                        tmp_celltype <- tail(split_name,1)
                        dataset_name <- paste0(paste(tmp_splitname,collapse=" ")," (",tmp_celltype,")")
                    }
                    Dataset <- c(Dataset, dataset_name)
                    NumSamples <- c(NumSamples, sample)
                    Perm <- c(Perm, str_split(perm,"_")[[1]][2])
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

    # create and return plot
    propDEGs_overall.plot <- ggplot(DEGs_df,aes(x=factor(numSamples),y=100*pctDEGs))+
                                    geom_boxplot(outlier.shape=NA,aes(fill=factor(Dataset)),position = position_dodge2(width = 1, preserve = "single"))+
                                    scale_y_continuous(labels = function(x) round(x))+
                                    guides(fill = guide_legend(label.format = function(x) round(as.numeric(x))))+
                                    theme_cowplot()+
                                    scale_fill_manual(values=generate_color_palette(length(list.dirs(inpath,recursive=F,full.names=F)),palette="Set1"))+
                                    labs(y="% DEGs", x="Number of Samples", fill="Dataset",title="Percentage DEGs from Bulk data detected when down-sampling scRNA-seq datasets (across all cell types)")+
                                    scale_alpha(guide = 'none')+
                                    theme(axis.title = element_text(size = fontsize_axislabels),
                                          axis.text.x = element_text(size = fontsize_axisticks),
                                          axis.text.y = element_text(size = fontsize_axisticks),
                                          plot.title = element_text(size = fontsize_title),
                                          legend.text = element_text(size = fontsize_legendlabels),
                                          legend.title = element_text(size = fontsize_legendtitle))

    return(propDEGs_overall.plot)

}