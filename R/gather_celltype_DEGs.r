#' Collate DEGs detected in DGE analysis outputs, across all celltypes in a dataset (datasets/DGE analysis outputs should have common celltypes as specified below)

#' @param celltype_correspondence list of different names specifying each cell type
#' @param pvalue the cut-off p-value used to select DEGs
#' @param Nperms number of permutations of DGE analysis outputs for each sample
#' @param output_path path storing the down-sampled DGE analysis for each single-cell dataset, generated for bulk analysis

#' Saves combined list of DEGs (across all cell types) in a subdirectory inside the dataset directory

gather_celltype_DEGs <- function(celltype_correspondence,
                                 pvalue=0.05,
                                 Nperms=20,
                                 output_path=getwd()){

    # validate function input params
    validate_input_parameters_bulk(output_path=output_path, celltype_correspondence=celltype_correspondence,
                                   pvalue=pvalue, Nperms=Nperms)

    # loop through datasets
    j <- 0
    for (dataset in list.dirs(output_path, recursive=FALSE, full.names=FALSE)) {
        j <- j + 1
        print(paste0("Gathering data for ", dataset))
        # base dataset path
        base_dataset <- file.path(output_path, dataset)
        # create path
        target_base <- file.path(base_dataset, "Overall/DE_downsampling")
        dir.create(target_base, recursive=TRUE, showWarnings=FALSE)
        # loop through each cell type
        for (standard_celltype in names(celltype_correspondence)) {
            # sampling point paths
            celltype_path <- file.path(base_dataset, standard_celltype, "DE_downsampling")
            downsampled_folders <- list.dirs(celltype_path, recursive=FALSE, full.names=FALSE)
            # loop through sampling points
            for (sample_samples in downsampled_folders) {
                # sample
                sample <- as.numeric(gsub("[^0-9]", "", sample_samples))
                print(paste0("Processing ", sample_samples, "..."))
                cell_sample_path <- file.path(celltype_path, sample_samples)
                # loop through permutations
                for (perm in 1:Nperms) {
                    # list to store DEGs
                    DEGs <- c()
                    sample_perm <- paste0(sample, "_", perm)
                    # create directory in overall folder
                    target_sample_perm <- file.path(target_base, sample_samples, sample_perm)
                    dir.create(target_sample_perm, recursive=TRUE, showWarnings=FALSE)
                    # check if the cell type directory exists
                    cell_sample_permpath <- file.path(cell_sample_path, sample_perm)
                    deg_file <- file.path(cell_sample_permpath, paste0("DEout", sample_perm, ".RData"))
                    if (file.exists(deg_file)) {
                        load(deg_file)
                        # get DEGs
                        DEGs <- c(DEGs, subset(eval(as.name(paste0("DEout_", sample)))$celltype_all_genes[[celltype_correspondence[[standard_celltype]][j]]], PValue <= pvalue)$name)
                        DEGs <- unique(DEGs)
                    }
                    # save DEGs
                    save(DEGs, file = file.path(target_sample_perm, paste0("DEout", sample_perm, ".RData")))
                }
            }
        }
    }
}
