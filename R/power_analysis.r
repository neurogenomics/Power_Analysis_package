# SCRIPT 4: Run Single DE Analysis and Stratify Results
#
# REFACTORED into a reusable function.
#
# GOAL:
# 1. Define a main function `run_de_analysis` that encapsulates
#    all the logic (load, run DE, stratify).
# 2. Call this function with a set of parameters to get a
#    "ground truth" DEG list.

# --- 1. SETUP ---
library(SingleCellExperiment)
library(qs)
library(scuttle)     # For aggregateAcrossCells
library(DESeq2)      # For the DE analysis
library(dplyr)       # For easy data filtering
library(tibble)      # For rownames_to_column

# --- 2. DEFINE THE REUSABLE DE ANALYSIS FUNCTION ---

#' Run Pseudobulk DE Analysis and Stratify Results
#'
#' Takes a single SCE object and performs a full DESeq2 analysis,
#' returning stratified lists of DEGs.
#'
#' @param sce A loaded SingleCellExperiment object.
#' @param sample_id_col The column name in colData for the sample/individual ID
#'   (e.g., "donor_id").
#' @param assay_name The string name of the assay to use (e.g., "counts" or "X").
#' @param design_formula A formula for the DESeq2 design (e.g., `~ sex + age`).
#'   All variables in the formula must exist in `colData(sce)`.
#' @param contrast_name Optional. The specific contrast to extract.
#'   If NULL (default), the function will attempt to auto-detect the
#'   primary contrast from the last variable in `design_formula`.
#'
#' @return A list containing three sub-lists:
#'   `user_method` (stratified by p-value and logFC),
#'   `standard_method` (stratified by adjusted p-value and logFC),
#'   `full_results` (the complete, non-stratified results data frame).
#'
run_de_analysis <- function(sce,
                            sample_id_col,
                            assay_name,
                            design_formula,
                            contrast_name = NULL) { # <-- PARAMETER IS NOW OPTIONAL

    message("Starting pseudobulk DE analysis...")

    # Check that the specified assay exists
    if (!assay_name %in% assayNames(sce)) {
        stop(paste0("Error: Assay '", assay_name, "' not found in SCE object. ",
                    "Available assays are: ", paste(assayNames(sce), collapse=", ")))
    }

    # 1. Aggregate counts to pseudobulk level
    message(paste("Aggregating by:", sample_id_col, "using assay:", assay_name))
    pb.sce <- aggregateAcrossCells(
        sce,
        ids = sce[[sample_id_col]],
        use.assay.type = assay_name # Pass the assay name here
    )

    # 2. Extract colData for DESeq2.
    pb.colData <- as.data.frame(colData(pb.sce))

    # --- Robust Variable Checking ---
    # Get all variables from the formula
    all_vars <- all.vars(design_formula)

    # Ensure all variables are in the colData and convert them to factors
    # This is safer than hard-coding 'sex'
    for (v in all_vars) {
        if (!v %in% colnames(pb.colData)) {
            stop(paste0("Error: Variable '", v, "' from design formula not in colData!"))
        }
        message(paste("Converting variable '", v, "' to factor for DE design."))
        pb.colData[[v]] <- as.factor(pb.colData[[v]])
    }
    # --- End Variable Checking ---

    # 3. Create the DESeqDataSet
    # *** Use assay_name parameter here ***
    dds <- DESeqDataSetFromMatrix(
        countData = assay(pb.sce, assay_name),
        colData = pb.colData,
        design = design_formula
    )

    # 4. Run the DESeq2 analysis
    message("Filtering low-count genes...")
    dds <- dds[rowSums(counts(dds)) >= 10, ] # Filter

    message("Running DESeq()...")
    dds <- DESeq(dds)

    # --- NEW: Automatic Contrast Detection ---

    # Get all available results names
    all_results_names <- resultsNames(dds)

    if (is.null(contrast_name)) {
        message("`contrast_name` not provided, attempting to auto-detect.")

        # Get the last variable from the design formula
        var_of_interest <- tail(all_vars, 1)

        # Find the first contrast that starts with this variable
        # This is the standard DESeq2 output for a factor
        found_contrasts <- grep(paste0("^", var_of_interest), all_results_names, value = TRUE)

        if (length(found_contrasts) == 0) {
            stop(paste0("Auto-detection failed. Could not find a contrast for variable '", var_of_interest, "'. ",
                        "Available contrasts are: ", paste(all_results_names, collapse=", "), ". ",
                        "Please specify one manually using the `contrast_name` parameter."))
        }

        # Use the first one found
        target_contrast <- found_contrasts[1]
        message(paste("Auto-detected contrast:", target_contrast))

    } else {
        message(paste("Using user-specified contrast:", contrast_name))
        target_contrast <- contrast_name
    }
    # --- END NEW SECTION ---


    # 5. Get the results
    message(paste("Extracting results for contrast:", target_contrast))

    if (!target_contrast %in% all_results_names) {
        warning(paste("Contrast '", target_contrast, "' not found. Available contrasts are: ",
                      paste(all_results_names, collapse=", ")))
        return(NULL)
    }

    res <- results(dds, name = target_contrast)
    res_df <- as.data.frame(res) %>%
        rownames_to_column("gene") %>%
        na.omit() # Remove NAs

    message("DE analysis complete. Stratifying results...")

    # --- 6. STRATIFY RESULTS (Using Adjusted P-Value as requested) ---
    # This section now implements the user's desired bins using 'padj'

    degs_stratified <- list()

    # Bin 1: padj < 0.05, |logFC| > 2
    degs_stratified$padj05_fc_gt_2 <- res_df %>%
        filter(padj < 0.05, abs(log2FoldChange) > 2) %>%
        arrange(padj)

    # Bin 2: padj < 0.05, 1 < |logFC| <= 2
    degs_stratified$padj05_fc_1_to_2 <- res_df %>%
        filter(padj < 0.05, abs(log2FoldChange) > 1, abs(log2FoldChange) <= 2) %>%
        arrange(padj)

    # Bin 3: padj < 0.05, 0 <= |logFC| <= 1
    degs_stratified$padj05_fc_0_to_1 <- res_df %>%
        filter(padj < 0.05, abs(log2FoldChange) >= 0, abs(log2FoldChange) <= 1) %>%
        arrange(padj)

    # Bin 4: "Buffer" list (using adjusted p-value)
    degs_stratified$padj_buffer_05_to_10 <- res_df %>%
        filter(padj >= 0.05, padj < 0.1) %>%
        arrange(padj)

    # Bin 5: All Significant (for a complete list, useful for power calculation)
    degs_stratified$all_significant_padj05 <- res_df %>%
        filter(padj < 0.05) %>%
        arrange(padj)

    # --- 8. RETURN ALL LISTS ---
    return(list(
        stratified_results = degs_stratified,
        full_results = res_df
    ))
}


# --- 3. EXAMPLE USAGE: GET "GROUND TRUTH" DEGS ---
#
# You can run this part of the script to get your baseline
# list of DEGs from the full, un-downsampled dataset.
#
message("\n--- Running 'Ground Truth' DE Analysis ---")

# --- Define Parameters ---
CELL_TYPE_DIR <- "/mnt/data/shared/poweranalysis/preprocessed_data_CONSOLIDATED/Astrocyte"
DATASET_FILE <- "Roussos_Combined.qs"
SAMPLE_COL <- "donor_id"
ASSAY_NAME <- "X" # *** ADDED THIS: Roussos data uses "X", not "counts" ***

# Define the DE design
# NOTE: Make sure the variables here (e.g., 'sex') exist in the colData!
FORMULA <- ~ sex
# CONTRAST <- "sex_M_vs_F" # <-- REMOVED to test auto-detection

# --- Load Ground Truth Data ---
message(paste("Loading:", file.path(CELL_TYPE_DIR, DATASET_FILE)))
ground_truth_sce <- qs::qread(file.path(CELL_TYPE_DIR, DATASET_FILE))

# --- Call the function ---
ground_truth_degs <- run_de_analysis(
    sce = ground_truth_sce,
    sample_id_col = SAMPLE_COL,
    assay_name = ASSAY_NAME,
    design_formula = FORMULA
    # contrast_name parameter is omitted, so auto-detection will run
)

# --- Print Summary of Ground Truth ---
if (!is.null(ground_truth_degs)) {
    message("\n--- Ground Truth Results (Adjusted P-Value Method) ---")
    message(paste("FDR < 0.05 & |logFC| > 2:",    nrow(ground_truth_degs$stratified_results$padj05_fc_gt_2)))
    message(paste("FDR < 0.05 & |logFC| 1-2:",  nrow(ground_truth_degs$stratified_results$padj05_fc_1_to_2)))
    message(paste("FDR < 0.05 & |logFC| 0-1:",  nrow(ground_truth_degs$stratified_results$padj05_fc_0_to_1)))
    message(paste("Buffer (0.05 < FDR < 0.1):",  nrow(ground_truth_degs$stratified_results$padj_buffer_05_to_10)))
    message(paste("Total Significant (FDR < 0.05):", nrow(ground_truth_degs$stratified_results$all_significant_padj05)))

    # You now have the 'ground_truth_degs' object.
    # You can save this object and use it in your main power analysis script
    # to compare against the results from downsampled datasets.

    # qs::qsave(ground_truth_degs, "ground_truth_astro_sex_degs.qs")
}

message("\n--- Script Complete ---")
