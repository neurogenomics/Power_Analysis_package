# Define global variables
utils::globalVariables(c("PValue","name"))

#' Downsample the dataset, based either on the individuals or cells, and run DE analysis on each downsampled output. Save results in a dataframe
#'
#' @importFrom stats as.formula
#'
#' @param SCE the input data (should be an SCE object)
#' @param range_downsampled vector or list containing values which the SCE will be downsampled at, in ascending order
#' @param output_path base path in which outputs will be stored
#' @param sampled downsampling carried out based on what (either "individuals" or "cells")
#' @param sampleID sample ID
#' @param design the design formula of class type `formula`. Equation used to fit the model- data for the generalised linear model e.g. expression ~ sex + pmi + disease
#' @param sexID sex ID
#' @param celltypeID cell type ID
#' @param assay_name the assay name in the SCE object to use for the analysis. Default is "counts" which uses the counts assay in each SCE
#' @param coef which coefficient to carry out DE analysis with respect to
#' @param fdr the cut-off False Discovery Rate below which to select DEGs
#' @param nom_pval the cut-off nominal P-value below which to select DEGs (as an alternative to FDR)
#' @param Nperms number of subsets created when downsampling at each level
#' @param y the column name in the SCE object for the return variable e.g. "diagnosis" - Case or disease. Default is the last variable in the design formula. y can be discrete (logistic regression) or continuous (linear regression)
#' @param region the column name in the SCE object for the study region. If there are multiple regions in the study (for example two brain regions). Pseudobulk values can be derived separately. Default is "single_region" which will not split by region.
#' @param control character specifying which control level for the differential expression analysis e.g. in a case/control/other study use "control" in the y column to compare against. NOTE only need to specify if more than two groups in y, leave as default value for two groups or continuous y. Default is NULL.
#' @param pval_adjust_method the adjustment method for the p-value in the differential expression analysis. Default is benjamini hochberg "BH". See  stats::p.adjust for available options
#' @param rmv_zero_count_genes whether genes with no count values in any cell should be removed. Default is TRUE
#'
#' @saves all DGE analysis outputs for downsampled files as well as a summary table of results showing number of true DEGs detected at each number of samples/cells

downsampling_DEanalysis <- function(SCE,
                                    range_downsampled="placeholder",
                                    output_path=getwd(),
                                    sampled="individuals",
                                    sampleID="donor_id",
                                    design="placeholder",
                                    sexID="sex",
                                    celltypeID="cell_type",
                                    assay_name="counts",
                                    coef="male",
                                    fdr=0.05,
                                    nom_pval=0.05,
                                    Nperms=20,
                                    y=NULL,
                                    region="single_region",
                                    control=NULL,
                                    pval_adjust_method="BH",
                                    rmv_zero_count_genes=TRUE){

    # --- 1. Initial Setup ---

    # alter range_downsampled
    if(identical(range_downsampled,"placeholder")){
        range_downsampled <- downsampling_range(SCE, sampled, sampleID)
    }
    # alter design
    if(design=="placeholder"){
        design=as.formula(paste0("~",sexID))
    }

    # create output path if doesn't already exist (no setwd())
    dir.create(output_path, showWarnings=FALSE, recursive=TRUE)

    # get celltype name from dataset
    celltype_name <- toString(unique(SCE[[celltypeID]]))

    # Define full path for ground truth file
    ground_truth_file <- file.path(output_path, "DEout.RData")

    # --- 2. Ground Truth Analysis ---
    # check if DE analysis output present already
    if(!file.exists(ground_truth_file)){
        message("Running full DGE analysis to establish ground truth...")
        # run and save DE analysis
        DEout <- DGE_analysis(SCE, design=design, sampleID=sampleID, celltypeID=celltypeID,
                              assay_name=assay_name, y=y, region=region, control=control,
                              pval_adjust_method=pval_adjust_method,
                              rmv_zero_count_genes=rmv_zero_count_genes,
                              verbose=T, coef=coef)
        save(DEout, file=ground_truth_file)
        message("Ground truth analysis complete.")
    } else {
        message("Loading existing ground truth DGE analysis from DEout.RData...")
        load(ground_truth_file) # This loads 'DEout' into the environment
    }

    # get DEGs using both FDR and nominal Pval
    all_genes_full <- DEout$celltype_all_genes[[celltype_name]]
    DEGs_fdr <- subset(all_genes_full, adj_pval < fdr)$name
    DEGs_pval <- subset(all_genes_full, PValue < nom_pval)$name

    # non DEGs
    nonDEGs_fdr <- subset(all_genes_full, !name %in% DEGs_fdr)$name
    nonDEGs_pval <- subset(all_genes_full, !name %in% nonDEGs_pval)$name

    # --- 3. Nested Helper Function ---
    # This function will run the inner loop (Nperms)
    # It is defined inside the main function so it can access all parameters
    # (design, coef, DEGs_fdr, etc.) without needing them as arguments.

    .run_perms_on_subsets <- function(subsets_list, value, path_val) {

        # Pre-allocate vectors to store results for this 'value'
        num_DEGs_fdr_vec <- numeric(Nperms)
        num_DEGs_pval_vec <- numeric(Nperms)
        num_FPs_fdr_vec <- numeric(Nperms)
        num_FPs_pval_vec <- numeric(Nperms)

        for(j in 1:Nperms){
            # create sub-directory for this permutation
            path_j <- file.path(path_val, paste0(toString(value),"_",j))
            dir.create(path_j, showWarnings=FALSE, recursive=TRUE)

            # Run DGE analysis on the j-th subset
            de_out_subset <- DGE_analysis(subsets_list[[j]], design=design, sampleID=sampleID,
                                          celltypeID=celltypeID, assay_name=assay_name, y=y,
                                          region=region, control=control,
                                          pval_adjust_method=pval_adjust_method,
                                          rmv_zero_count_genes=rmv_zero_count_genes,
                                          verbose=T, coef=coef)

            # save output
            save(de_out_subset, file=file.path(path_j, paste0("DEout",toString(value),"_",j,".RData")))

            # get number of TP DEGs
            all_genes_subset <- de_out_subset$celltype_all_genes[[celltype_name]]
            degs_fdr <- subset(all_genes_subset, adj_pval < fdr)$name
            degs_pval <- subset(all_genes_subset, PValue < nom_pval)$name

            num_DEGs_fdr_vec[j] <- sum(degs_fdr %in% DEGs_fdr)
            num_DEGs_pval_vec[j] <- sum(degs_pval %in% DEGs_pval)

            # get number of FP DEGs
            num_FPs_fdr_vec[j] <- sum(degs_fdr %in% nonDEGs_fdr)
            num_FPs_pval_vec[j] <- sum(degs_pval %in% nonDEGs_pval)
        }

        # Return a list of the completed vectors
        return(list(
            tp_fdr = num_DEGs_fdr_vec,
            tp_pval = num_DEGs_pval_vec,
            fp_fdr = num_FPs_fdr_vec,
            fp_pval = num_FPs_pval_vec
        ))
    }

    # --- 4. Main Downsampling Loop ---

    # Initialize lists to store result vectors for all 'values'
    all_results <- list(
        DEGs_detected_fdr = list(),
        DEGs_detected_pval = list(),
        falsePositives_fdr = list(),
        falsePositives_pval = list()
    )

    final_output_path <- "" # Will be set inside the loop

    for(value in range_downsampled){

        message(paste("--- Starting downsampling for value:", value, "---"))

        # 1. Set up variables and paths based on 'sampled' type
        if(sampled == "individuals"){
            output_dir_name <- "DE_downsampling"
            dir_suffix <- "samples"
            col_suffix <- "samples"

            # Create subsets
            subsets <- sample_individuals(SCE, value, sampleID, sexID, Nperms)

        } else {
            output_dir_name <- "DE_downsampling_cells"
            dir_suffix <- "cells_persample"
            col_suffix <- "cells_persample"

            # Create subsets (note the different arguments from sample_individuals)
            subsets <- sample_cells(SCE, value, sampleID, Nperms)
        }

        # Set and create paths
        final_output_path <- file.path(output_path, output_dir_name)
        path_val <- file.path(final_output_path, paste0(toString(value), dir_suffix))
        dir.create(path_val, showWarnings=FALSE, recursive=TRUE)

        # 2. Run helper function on the list of subsets
        perm_results <- .run_perms_on_subsets(subsets, value, path_val)

        # 3. Store results in the main lists using named elements
        tp_fdr_col_name <- paste0("numDEGs_", toString(value), col_suffix)
        tp_pval_col_name <- paste0("numDEGs_", toString(value), col_suffix)
        fp_fdr_col_name <- paste0("numFPs_", toString(value), col_suffix)
        fp_pval_col_name <- paste0("numFPs_", toString(value), col_suffix)

        all_results$DEGs_detected_fdr[[tp_fdr_col_name]] <- perm_results$tp_fdr
        all_results$DEGs_detected_pval[[tp_pval_col_name]] <- perm_results$tp_pval
        all_results$falsePositives_fdr[[fp_fdr_col_name]] <- perm_results$fp_fdr
        all_results$falsePositives_pval[[fp_pval_col_name]] <- perm_results$fp_pval

        message(paste("--- Finished for value:", value, "---"))
    } # End of main 'for' loop


    # --- 5. Save Final Results ---

    message("Downsampling complete. Saving summary data frames...")

    # Convert lists of vectors into data frames
    Iteration <- 1:Nperms
    DEGs_detected_fdr <- as.data.frame(all_results$DEGs_detected_fdr)
    DEGs_detected_fdr$Iteration <- Iteration

    DEGs_detected_pval <- as.data.frame(all_results$DEGs_detected_pval)
    DEGs_detected_pval$Iteration <- Iteration

    falsePositives_fdr <- as.data.frame(all_results$falsePositives_fdr)
    falsePositives_fdr$Iteration <- Iteration

    falsePositives_pval <- as.data.frame(all_results$falsePositives_pval)
    falsePositives_pval$Iteration <- Iteration

    # Save the final data frames to the correct sub-directory
    save(DEGs_detected_fdr, file = file.path(final_output_path, "DEGs_detected_fdr.RData"))
    save(DEGs_detected_pval, file = file.path(final_output_path, "DEGs_detected_pval.RData"))
    save(falsePositives_fdr, file = file.path(final_output_path, "falsePositives_fdr.RData"))
    save(falsePositives_pval, file = file.path(final_output_path, "falsePositives_pval.RData"))

    message(paste("All results saved to:", final_output_path))
}
