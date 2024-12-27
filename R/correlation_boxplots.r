# Define global variables
utils::globalVariables(c("alpha"))

#' Obtain box plots for the correlations of all celltypes, and the mean correlations at a specified cutoff p-value

#' @importFrom ggplot2 geom_boxplot ggplot labs facet_wrap theme aes element_text element_blank scale_colour_brewer scale_alpha guides guide_legend
#' @importFrom cowplot theme_cowplot

#' @param corrMats (named) list of correlation matrices for each celltype with the final element being the mean correlation matrix, all at specified p-value
#' @param numRealDatasets total number of *real* datasets (most likely the number of studies, but sometimes a study may be split e.g. into 2 brain regions, so in this case it would be the number of studies plus 1)
#' @param pval the cut-off p-value which was used to select DEGs
#' @param alphaval (alpha) transparency of the non-mean boxplots
#' @param numPerms number of random permutations of the dataset used to select significant DEGs from
#' @param numSubsets number of pairs of random subsets of the dataset used to select significant DEGs from
#' @param sexDEGs true if DEGs come from sex chromosomes, else false
#' @param fontsize_yaxislabels font size for axis labels in plot
#' @param fontsize_yaxisticks font size for axis tick labels in plot
#' @param fontsize_title font size for plot title
#' @param fontsize_legendlabels font size for legend labels in plot
#' @param fontsize_legendtitle font size for legend title in plot
#' @param fontsize_facet_labels font size for facet labels

#' @return box plots for correlation matrices at a certain p-value cut-off, sorted by celltype and then type of correlation

correlation_boxplots <- function(corrMats,
                                 numRealDatasets,
                                 pval,
                                 alphaval=0.25,
                                 numPerms=5,
                                 numSubsets=5,
                                 sexDEGs=FALSE,
                                 fontsize_yaxislabels=12,
                                 fontsize_yaxisticks=9,
                                 fontsize_title=14,
                                 fontsize_legendlabels=9,
                                 fontsize_legendtitle=9,
                                 fontsize_facet_labels=9){

    # validate function input params
    validate_input_parameters_correlation(corrMats=corrMats, numRealDatasets=numRealDatasets, pvalue=pval,
                                          alphaval=alphaval, numPerms=numPerms, numSubsets=numSubsets,
                                          sexDEGs=sexDEGs, fontsize_yaxislabels=fontsize_yaxislabels, fontsize_yaxisticks=fontsize_yaxisticks,
                                          fontsize_title=fontsize_title, fontsize_legendlabels=fontsize_legendlabels, fontsize_legendtitle=fontsize_legendtitle,
                                          fontsize_facet_labels=fontsize_facet_labels)

    # midCor submatrix limits
    midCorLim <- numPerms + numRealDatasets
    # list to hold results
    corrOuts <- c()
    # index
    j <- 1

    # get lists with all correlations
    for(corrMat in corrMats){
        # specify submatrices with upper/middle/lower bounds
        lower <- corrMat[1:numPerms,1:numPerms]
        middle <- corrMat[(numPerms+1):midCorLim,(numPerms+1):midCorLim]
        upper <- corrMat[(midCorLim+1):(midCorLim+numSubsets),(midCorLim+1):(midCorLim+numSubsets)]
        # convert each one to a list and remove "1" (selfcorrelation)
        lower <- unique(unlist(as.list(lower)))
        lower <- lower[-c(1)]
        middle <- unique(unlist(as.list(middle)))
        middle <- middle[-c(1)]
        upper <- unique(unlist(as.list(upper)))
        upper <- upper[-c(1)]
        # store in list
        corrOuts[[j]] <- list(lower,middle,upper)
        names(corrOuts)[[j]] <- names(corrMats)[[j]]
        # increment
        j <- j+1
    }

    # store in dataframe
    i <- 1
    # empty dataframe
    df <- data.frame()
    # fill dataframe
    for(out in corrOuts){
        # define variables
        var1 <- replicate(length(out[[1]])+length(out[[2]])+length(out[[3]]), names(corrOuts)[[i]]) #celltype
        var2 <- c(replicate(length(out[[1]]),"Random Permutations"),replicate(length(out[[2]]),"Between Study"),replicate(length(out[[3]]),"Within-study subsamples"))
        val <- unlist(out)
        # put in dataframe
        df_new <- data.frame(var1)
        df_new$var2 <- var2
        df_new$val <- val
        # join
        df <- rbind(df,df_new)
        i <- i+1
    }

    df$alpha <- replicate(dim(df)[[1]],alphaval)
    df$alpha <- ifelse(df$var1=="Mean", 1, ifelse(df$var1!="Mean",alphaval,alphaval))
    unique_alphas <- df[!duplicated(df[,c("var1")]),]$alpha

    # box plot
    if(pval == 1){
        if(sexDEGs == FALSE){
            fig.plot <- ggplot(df,
                aes(x=factor(var1,levels=c("Mean","Astro","Endo","Micro","Oligo")),y=val))+
                geom_boxplot(outlier.shape=NA,aes(fill=factor(var1,levels=c("Mean","Astro","Endo","Micro","Oligo")),alpha=alpha),width=5)+
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
                aes(x=factor(var1,levels=c("Mean","Astro","Endo","Micro","Oligo")),y=val))+
                geom_boxplot(outlier.shape=NA,aes(fill=factor(var1,levels=c("Mean","Astro","Endo","Micro","Oligo")),alpha=alpha),width=5)+
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
        if(sexDEGs == FALSE){
            fig.plot <- ggplot(df,
                aes(x=factor(var1,levels=c("Mean","Astro","Endo","Micro","Oligo")),y=val))+
                geom_boxplot(outlier.shape=NA,aes(fill=factor(var1,levels=c("Mean","Astro","Endo","Micro","Oligo")),alpha=alpha),width=5)+
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
                aes(x=factor(var1,levels=c("Mean","Astro","Endo","Micro","Oligo")),y=val))+
                geom_boxplot(outlier.shape=NA,aes(fill=factor(var1,levels=c("Mean","Astro","Endo","Micro","Oligo")),alpha=alpha),width=5)+
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

    return(fig.plot)

}