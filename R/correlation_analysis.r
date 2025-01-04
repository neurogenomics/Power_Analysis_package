#' Runs correlation analysis pipeline

#' @param dataset_name name of the dataset used to select significant DEGs from (specified as a string, name as in allstudies)
#' @param allstudies a list containing all the datasets (as SCE objects)
#' @param celltypes a list containing the celltypes to compute mean correlation across
#' @param pvals list of p-value cut-offs which will be used to select DEGs
#' @param data_names names of the datasets as they appear in the correlation plot
#' @param alphaval (alpha) transparency of the non-mean boxplots
#' @param numPerms number of random permutations of the dataset used to select significant DEGs from
#' @param numSubsets number of pairs of random subsets of the dataset used to select significant DEGs from
#' @param output_path base path in which outputs will be stored

#' Saves all plots and DE outputs in the appropriate directories
#' @export

correlation_analysis <- function(dataset_name="placeholder",
                                 allstudies="placeholder",
                                 celltypes="placeholder",
                                 pvals=c(0.05,0.025,0.01,0.001,0.0001),
                                 data_names="placeholder",
                                 alphaval="placeholder",
                                 numPerms="placeholder",
                                 numSubsets="placeholder",
                                 output_path=getwd()){
    
    # run plot_mean_correlation for each p-value (saving outputs)
    mean_correlation_results <- plot_mean_correlation(dataset_name,
                                                      allstudies,
                                                      celltypes,
                                                      pvals,
                                                      data_names,
                                                      output_path)

    # run correlation_boxplots for each p-value (saving outputs)
    boxplot_results <- correlation_boxplots(mean_correlation_results,
                                            numRealDatasets=length(allstudies),
                                            pvals=pvals,
                                            output_path=output_path)

    return(list(mean_correlation_results = mean_correlation_results, boxplot_results = boxplot_results))

}