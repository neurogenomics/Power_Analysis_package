#' Collate DEGs detected in DGE analysis outputs, across all celltypes in a dataset (datasets/DGE analysis outputs should have common celltypes as specified below)

#' @param celltype_correspondence A named vector that maps a standard cell type label (e.g., `"Endo"`, `"Micro"`) to how that cell type appears in each dataset. Use `NA` if the cell type is not present in a given dataset.
#' @param pvalue P-value threshold for defining DEGs in the bulk dataset.
#' @param Nperms Number of permutations to perform for each down-sampling level. Default is 20.
#' @param output_path A clean directory path where down-sampled outputs and plots will be saved (should contain no subdirectories).
#' @param sampled Specifies the unit of down-sampling. Can be either `"individuals"` or `"cells"`, depending on whether the analysis downsamples across samples or cells.

#' Saves combined list of DEGs (across all cell types) in a subdirectory inside the dataset directory

gather_celltype_DEGs <- function(celltype_correspondence,
                                 pvalue=0.05,
                                 Nperms=20,
                                 sampled="individuals",
                                 output_path=getwd()){

    sampled <- match.arg(sampled, choices = c("individuals", "cells"))

    # validate function input params
    validate_input_parameters_bulk(output_path=output_path, celltype_correspondence=celltype_correspondence,
                                   pvalue=pvalue, Nperms=Nperms)

    # reference the appropriate subdirectory regarding the down-sampling type
    folder_tag <- if (sampled == "individuals") "DE_downsampling" else "DE_downsampling_cells"

    # loop through datasets
    j <- 0
    for (dataset in list.dirs(output_path, recursive=FALSE, full.names=FALSE)) {
        j <- j + 1
        print(paste0("Gathering data for ", dataset))
        # base dataset path
        base_dataset <- file.path(output_path, dataset)
        # create path
        target_base <- file.path(base_dataset, "Overall",folder_tag)
        dir.create(target_base, recursive=TRUE, showWarnings=FALSE)
        # loop through each cell type
        for (standard_celltype in names(celltype_correspondence)) {
            # sampling point paths
            celltype_path <- file.path(base_dataset, standard_celltype, folder_tag)
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
