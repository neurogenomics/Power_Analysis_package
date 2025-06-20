# Define global variables
utils::globalVariables(c("alpha"))

#' Obtain box plots for the correlations of all celltypes, and the mean correlations at a specified cutoff p-value

#' @importFrom ggplot2 geom_boxplot ggplot labs facet_wrap theme aes element_text element_blank scale_colour_brewer scale_alpha guides guide_legend
#' @importFrom cowplot theme_cowplot

#' @param corr_mats (named) list of correlation matrices for each celltype with the final element being the mean correlation matrix, all at specified p-value
#' @param num_real_datasets total number of *real* datasets (most likely the number of studies, but sometimes a study may be split e.g. into 2 brain regions, so in this case it would be the number of studies plus 1)
#' @param pvals the cut-off p-value which was used to select DEGs
#' @param alphaval (alpha) transparency of the non-mean boxplots
#' @param N_randperms number of random permutations of the dataset used to select significant DEGs from
#' @param N_subsets number of pairs of random subsets of the dataset used to select significant DEGs from
#' @param sex_DEGs true if DEGs come from sex chromosomes, else false
#' @param fontsize_yaxislabels font size for axis labels in plot
#' @param fontsize_yaxisticks font size for axis tick labels in plot
#' @param fontsize_title font size for plot title
#' @param fontsize_legendlabels font size for legend labels in plot
#' @param fontsize_legendtitle font size for legend title in plot
#' @param fontsize_facet_labels font size for facet labels
#' @param output_path base path in which outputs will be stored

#' Saves box plots for correlation matrices, at a certain p-value cut-off, in the appropriate directory 

correlation_boxplots <- function(corr_mats,
                                 num_real_datasets,
                                 pvals=c(0.05,0.025,0.01,0.001,0.0001),
                                 alphaval=0.25,
                                 N_randperms=5,
                                 N_subsets=5,
                                 sex_DEGs=FALSE,
                                 fontsize_yaxislabels=12,
                                 fontsize_yaxisticks=9,
                                 fontsize_title=14,
                                 fontsize_legendlabels=9,
                                 fontsize_legendtitle=9,
                                 fontsize_facet_labels=9,
                                 output_path=getwd()){                                    
    # outputs
    output_plots <- list()

    # loop over each p-value
    for(pval in pvals){
        # midCor submatrix limits
        midCorLim <- N_randperms + num_real_datasets
        # correlation matrix for this p-value
        corr_mats_pval <- corr_mats[[as.character(pval)]]
        # df to hold results
        df <- data.frame()
        # index
        j <- 1

        # get lists with all correlations
        for(celltype in names(corr_mats_pval)){
            # skip if invalid matrix
            corrMat <- corr_mats_pval[[celltype]]
            # check if corrMat is a valid matrix
            if (!is.matrix(corrMat) || length(dim(corrMat)) != 2) {
                warning(paste0("Invalid correlation matrix for cell type ", celltype, " at p-value ", pval, ". Skipping."))
                next
            }
            # specify submatrices with upper/middle/lower bounds
            lower <- corrMat[1:N_randperms,1:N_randperms]
            middle <- corrMat[(N_randperms+1):midCorLim,(N_randperms+1):midCorLim]
            upper <- corrMat[(midCorLim+1):(midCorLim+N_subsets),(midCorLim+1):(midCorLim+N_subsets)]
            # convert each one to a list and remove "1" (self-correlation)
            lower <- unique(unlist(as.list(lower)))
            lower <- lower[-c(1)]
            middle <- unique(unlist(as.list(middle)))
            middle <- middle[-c(1)]
            upper <- unique(unlist(as.list(upper)))
            upper <- upper[-c(1)]
            # variables for df
            var1 <- replicate(length(lower) + length(middle) + length(upper), celltype)
            var2 <- c(replicate(length(lower), "Random Permutations"), replicate(length(middle), "Between Study"), replicate(length(upper), "Within-study subsamples"))
            val <- c(lower, middle, upper)
            # create a new dataframe for the current cell type
            df_new <- data.frame(var1 = var1, var2 = var2, val = val)
            # append to the main dataframe
            df <- rbind(df, df_new)
        }
        # add alpha col.
        df$alpha <- replicate(dim(df)[[1]],alphaval)
        df$alpha <- ifelse(df$var1=="Mean", 1, alphaval)
        unique_alphas <- df[!duplicated(df[,c("var1")]),]$alpha

        # check if df is valid (non-empty and contains required columns)
        if (nrow(df) == 0 || !all(c("var1", "var2", "val") %in% colnames(df))) {
            print(paste0("Data frame 'df' is invalid or empty for p-value ", pval, ". Skipping plot creation."))
            next
        }
        # levels
        unique_levels <- unique(df$var1)
        unique_levels <- c("Mean", setdiff(unique_levels, "Mean"))  # ensure "Mean" is always first
        # title
        plt_title <- ifelse(pval==1 && sex_DEGs==FALSE, "Using all DEGs (from all chromosomes)",
                     ifelse(pval==1 && sex_DEGs==TRUE, "Using all DEGs (from sex chromosomes)",
                     ifelse(pval!=1 && sex_DEGs==FALSE, paste0("DEGs selected at a ", pval*100, "% cut-off (from all chromosomes)"),
                     ifelse(pval!=1 && sex_DEGs==TRUE, paste0("DEGs selected at a ", pval*100, "% cut-off (from sex chromosomes)")))))
        # box plot
        fig.plot <- ggplot(df,
            aes(x=factor(var1,levels=unique_levels),y=val))+
            geom_boxplot(outlier.shape=NA,aes(fill=factor(var1,levels=unique_levels),alpha=alpha),width=5)+
            theme_cowplot()+
            scale_colour_brewer(palette = "Set1")+
            labs(y="Correlation", x="Type", fill="Cell Type",title=plt_title)+
            facet_wrap(factor(df$var2,levels=c("Random Permutations","Between Study","Within-study subsamples")), scales="fixed")+
            theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.title.y = element_text(size = fontsize_yaxislabels),
                axis.text.y = element_text(size = fontsize_yaxisticks),
                plot.title = element_text(size = fontsize_title),
                legend.text = element_text(size = fontsize_legendlabels),
                legend.title = element_text(size = fontsize_legendtitle),
                strip.text = element_text(size = fontsize_facet_labels))+
            scale_alpha(guide = 'none')+
            guides(fill=guide_legend(override.aes = list(alpha=unique_alphas)))

        # store output in list with p-value as key
        output_plots[[as.character(pval)]] <- fig.plot

        # save the plot if output_path is specified
        ggsave(file.path(output_path,paste0("correlation_boxplot_p", pval, ".png")), fig.plot, width=30, height=15, units="cm", bg="white")
        ggsave(file.path(output_path,paste0("correlation_boxplot_p", pval, ".pdf")), fig.plot, width=30, height=15, units="cm", bg="white")

    }

}