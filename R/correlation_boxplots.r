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

    # validate function input params
    validate_input_parameters_correlation(corr_mats=corr_mats, num_real_datasets=num_real_datasets, pvalues=pvals,
                                          alphaval=alphaval, N_randperms=N_randperms, N_subsets=N_subsets,
                                          sex_DEGs=sex_DEGs, fontsize_yaxislabels=fontsize_yaxislabels, fontsize_yaxisticks=fontsize_yaxisticks,
                                          fontsize_title=fontsize_title, fontsize_legendlabels=fontsize_legendlabels, fontsize_legendtitle=fontsize_legendtitle,
                                          fontsize_facet_labels=fontsize_facet_labels, output_path=output_path)                                    
    # outputs
    output_plots <- list()

    # loop over each p-value
    for(pval in pvals){
        # midCor submatrix limits
        midCorLim <- N_randperms + num_real_datasets
        # list to hold results
        corrOuts <- c()
        # index
        j <- 1

        # get lists with all correlations
        if(length(corr_mats) > 0){
            for(corrMat in corr_mats){
                if (!is.matrix(corrMat) || length(dim(corrMat)) != 2) {
                    print("Invalid correlation matrix encountered - skipping.")
                    next
                }
                # specify submatrices with upper/middle/lower bounds
                lower <- corrMat[1:N_randperms,1:N_randperms]
                middle <- corrMat[(N_randperms+1):midCorLim,(N_randperms+1):midCorLim]
                upper <- corrMat[(midCorLim+1):(midCorLim+N_subsets),(midCorLim+1):(midCorLim+N_subsets)]
                # convert each one to a list and remove "1" (selfcorrelation)
                lower <- unique(unlist(as.list(lower)))
                lower <- lower[-c(1)]
                middle <- unique(unlist(as.list(middle)))
                middle <- middle[-c(1)]
                upper <- unique(unlist(as.list(upper)))
                upper <- upper[-c(1)]
                # store in list
                corrOuts[[j]] <- list(lower,middle,upper)
                names(corrOuts)[[j]] <- names(corr_mats)[[j]]
                # increment
                j <- j+1
            }
        }else{
            # if no correlation matrices are provided, create empty lists
            corrOuts <- list()
        }
        # store in dataframe
        i <- 1
        # empty dataframe
        df <- data.frame()
        if(length(corrOuts > 0)){
            # fill dataframe
            for(out in corrOuts){
                # define variables
                var1 <- replicate(length(out[[1]])+length(out[[2]])+length(out[[3]]), names(corrOuts)[[i]]) #celltype
                var2 <- c(replicate(length(out[[1]]),"Random Permutations"),replicate(length(out[[2]]),"Between Study"),replicate(length(out[[3]]),"Within-study subsamples"))
                val <- unlist(out)
                # put in dataframe
                df_new <- data.frame(var1=var1, var2=var2, val=val)
                # join
                df <- rbind(df,df_new)
                i <- i+1
            }
        }else{
            # if no correlation matrices are provided, create empty dataframe
            df <- data.frame(var1=character(), var2=character(), val=numeric())
        }
        df$alpha <- replicate(dim(df)[[1]],alphaval)
        df$alpha <- ifelse(df$var1=="Mean", 1, ifelse(df$var1!="Mean",alphaval,alphaval))
        unique_alphas <- df[!duplicated(df[,c("var1")]),]$alpha

        # check if df is valid (non-empty and contains required columns)
        if (nrow(df) == 0 || !all(c("var1", "var2", "val") %in% colnames(df))) {
            print(paste0("Data frame 'df' is invalid or empty for p-value ", pval, ". Skipping plot creation."))
            next
        }
        # levels
        unique_levels <- unique(df$var1)
        unique_levels <- c("Mean", setdiff(unique_levels, "Mean"))  # ensure "Mean" is always first
        # box plot
        if(pval == 1){
            if(sex_DEGs == FALSE){
                fig.plot <- ggplot(df,
                    aes(x=factor(var1,levels=unique_levels),y=val))+
                    geom_boxplot(outlier.shape=NA,aes(fill=factor(var1,levels=unique_levels),alpha=alpha),width=5)+
                    theme_cowplot()+
                    scale_colour_brewer(palette = "Set1")+
                    labs(y="Correlation", x="Type", fill="Cell Type",title=paste0("Using all DEGs (from all chromosomes)"))+
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
            }else{
                fig.plot <- ggplot(df,
                    aes(x=factor(var1,levels=unique_levels),y=val))+
                    geom_boxplot(outlier.shape=NA,aes(fill=factor(var1,levels=unique_levels),alpha=alpha),width=5)+
                    theme_cowplot()+
                    scale_colour_brewer(palette = "Set1")+
                    labs(y="Correlation", x="Type", fill="Cell Type",title=paste0("Using all DEGs (from sex chromosomes)"))+
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
            }
        }else{
            if(sex_DEGs == FALSE){
                fig.plot <- ggplot(df,
                    aes(x=factor(var1,levels=unique_levels),y=val))+
                    geom_boxplot(outlier.shape=NA,aes(fill=factor(var1,levels=unique_levels),alpha=alpha),width=5)+
                    theme_cowplot()+
                    scale_colour_brewer(palette = "Set1")+
                    labs(y="Correlation", x="Type", fill="Cell Type",title=paste0("DEGs selected at a ",pval*100,"% cut-off (from all chromosomes)"))+
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
            }else{
                fig.plot <- ggplot(df,
                    aes(x=factor(var1,levels=unique_levels),y=val))+
                    geom_boxplot(outlier.shape=NA,aes(fill=factor(var1,levels=unique_levels),alpha=alpha),width=5)+
                    theme_cowplot()+
                    scale_colour_brewer(palette = "Set1")+
                    labs(y="Correlation", x="Type", fill="Cell Type",title=paste0("DEGs selected at a ",pval*100,"% cut-off (from sex chromosomes)"))+
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
            }
        }

        # store output in list with p-value as key
        output_plots[[as.character(pval)]] <- fig.plot

        # save the plot if output_path is specified
        ggsave(file.path(output_path,paste0("correlation_boxplot_p", pval, ".png")), fig.plot, width=30, height=15, units="cm", bg="white")
        ggsave(file.path(output_path,paste0("correlation_boxplot_p", pval, ".pdf")), fig.plot, width=30, height=15, units="cm", bg="white")

    }

}