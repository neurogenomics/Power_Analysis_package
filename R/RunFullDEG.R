# SCRIPT 4: Run Single DE Analysis (PARALLEL + GRANULAR CHECKPOINTS + FAST MODE)
#
# GOAL:
# 1. Define a main function `run_de_analysis` that encapsulates
#    all the logic (load, run DE, stratify).
# 2. Call this function with a set of parameters to get a
#    "ground truth" DEG list.

# --- 1. SETUP ---
library(SingleCellExperiment)
library(qs)
library(scuttle)       # For aggregateAcrossCells
library(DESeq2)        # For the DE analysis
library(dplyr)         # For easy data filtering
library(tibble)        # For rownames_to_column
library(BiocParallel)  # *** ADDED for parallelization ***

# --- REGISTER PARALLEL BACKEND (16 Cores) ---
N_CORES <- 32
message(paste("Registering BiocParallel backend with", N_CORES, "cores."))
register(MulticoreParam(workers = N_CORES))
# --------------------------------------------------


# --- 2. DEFINE THE REUSABLE DE ANALYSIS FUNCTION ---

#' Run Pseudobulk DE Analysis and Stratify Results
#'
#' Takes a single SCE object and performs a full DESeq2 analysis,
#' returning stratified lists of DEGs.
#'
#' @param sce A loaded SingleCellExperiment object.
#' @param sample_id_col Column for sample/individual ID.
#' @param assay_name String name of the assay to use.
#' @param design_formula A formula for the DESeq2 design.
#' @param contrast_name Optional. The specific contrast to extract.
#' @param pb_checkpoint Optional. Path to save/load the aggregated pseudobulk object.
#' @param dds_disp_checkpoint Optional. Path for dispersion checkpoint.
#' @param dds_wald_checkpoint Optional. Path for final fitted (Wald test) checkpoint.
#' @param fast_subset_cells Optional. If set to a number (e.g., 50000),
#'     subsamples the SCE to this many cells for a fast test run.
#'
#' @return A list containing stratified DEGs and full results.
#'
run_de_analysis <- function(sce,
                            sample_id_col,
                            assay_name,
                            design_formula,
                            contrast_name = NULL,
                            pb_checkpoint = NULL,
                            dds_disp_checkpoint = NULL,
                            dds_wald_checkpoint = NULL,
                            fast_subset_cells = NULL) { # <-- NEW PARAMETER

    message("Starting pseudobulk DE analysis...")

    # --- *** NEW: FAST MODE SUBSETTING *** ---
    if (!is.null(fast_subset_cells) && fast_subset_cells < ncol(sce)) {
        message(paste("--- FAST MODE: Subsampling to", fast_subset_cells, "random cells ---"))
        # Take a random sample of cells
        sce <- sce[, sample(ncol(sce), fast_subset_cells)]
    }
    # -------------------------------------------

    # --- *** 1. AGGREGATION CHECKPOINT *** ---
    if (!is.null(pb_checkpoint) && file.exists(pb_checkpoint)) {
        message(paste("Loading aggregated data from checkpoint:", pb_checkpoint))
        pb.sce <- qs::qread(pb_checkpoint)
    } else {
        # ... (Aggregation logic) ...
        message("Checkpoint not found. Running aggregation...")
        if (!assay_name %in% assayNames(sce)) {
            stop(paste0("Error: Assay '", assay_name, "' not found."))
        }
        if (any(is.na(sce[[sample_id_col]]))) {
            n_na <- sum(is.na(sce[[sample_id_col]]))
            message(paste("Warning: Removing", n_na, "cells with NA donor_id."))
            sce <- sce[, !is.na(sce[[sample_id_col]])]
        }
        message(paste("Aggregating", ncol(sce), "cells by:", sample_id_col, "..."))
        pb.sce <- aggregateAcrossCells(
            sce, ids = sce[[sample_id_col]], use.assay.type = assay_name
        )
        message(paste("Aggregation complete. Dimensions:",
                      nrow(pb.sce), "genes x", ncol(pb.sce), "individuals."))
        if (!is.null(pb_checkpoint)) {
            message(paste("Saving aggregation checkpoint to:", pb_checkpoint))
            qs::qsave(pb.sce, pb_checkpoint)
        }
    }

    # --- *** CRITICAL MEMORY MANAGEMENT *** ---
    if (exists("sce", inherits = FALSE)) {
        message("Removing large original SCE object from memory...")
        rm(sce)
        gc()
        message("Memory freed. Proceeding to DESeq2.")
    }

    # --- *** 2. DESEQ2 ANALYSIS (GRANULAR CHECKPOINTS) *** ---

    if (!is.null(dds_wald_checkpoint) && file.exists(dds_wald_checkpoint)) {
        message(paste("Loading *final* fitted DESeq object from checkpoint:", dds_wald_checkpoint))
        dds <- qs::qread(dds_wald_checkpoint)
    } else {
        if (!is.null(dds_disp_checkpoint) && file.exists(dds_disp_checkpoint)) {
            message(paste("Loading *dispersion* fitted object from checkpoint:", dds_disp_checkpoint))
            dds <- qs::qread(dds_disp_checkpoint)
        } else {
            message("No DESeq checkpoints found. Building DESeqDataSet...")
            pb.colData <- as.data.frame(colData(pb.sce))
            dds <- DESeqDataSetFromMatrix(
                countData = assay(pb.sce, assay_name),
                colData = pb.colData,
                design = design_formula
            )
            message(paste("Filtering low-count genes (initial count:", nrow(dds), ")..."))
            dds <- dds[rowSums(counts(dds)) >= 10, ]
            message(paste("Genes remaining after filtering:", nrow(dds)))

            message("Estimating size factors...")
            dds <- estimateSizeFactors(dds)

            message(paste("Estimating dispersions in parallel (", N_CORES, " cores)..."))
            dds <- estimateDispersions(dds, parallel = TRUE)
            message("Dispersion estimation complete.")

            if (!is.null(dds_disp_checkpoint)) {
                message(paste("Saving dispersion checkpoint to:", dds_disp_checkpoint))
                qs::qsave(dds, dds_disp_checkpoint)
            }
        }

        message(paste("Running Wald test in parallel (", N_CORES, " cores)..."))
        dds <- nbinomWaldTest(dds, parallel = TRUE)
        message("Wald test complete.")

        if (!is.null(dds_wald_checkpoint)) {
            message(paste("Saving final fitted (Wald) checkpoint to:", dds_wald_checkpoint))
            qs::qsave(dds, dds_wald_checkpoint)
        }
    }

    # --- (Fast steps: Contrast detection and result stratification) ---
    all_results_names <- resultsNames(dds)
    all_vars <- all.vars(design_formula)
    if (is.null(contrast_name)) {
        var_of_interest <- tail(all_vars, 1)
        found_contrasts <- grep(paste0("^", var_of_interest), all_results_names, value = TRUE)
        target_contrast <- found_contrasts[1]
    } else {
        target_contrast <- contrast_name
    }
    message(paste("Extracting results for contrast:", target_contrast))
    res <- results(dds, name = target_contrast)
    res_df <- as.data.frame(res) %>% rownames_to_column("gene") %>% na.omit()

    # ... (Stratification logic remains the same) ...
    degs_stratified <- list()
    degs_stratified$padj05_fc_gt_2 <- res_df %>%
        filter(padj < 0.05, abs(log2FoldChange) > 2) %>% arrange(padj)
    degs_stratified$padj05_fc_1_to_2 <- res_df %>%
        filter(padj < 0.05, abs(log2FoldChange) > 1, abs(log2FoldChange) <= 2) %>% arrange(padj)
    degs_stratified$padj05_fc_0_to_1 <- res_df %>%
        filter(padj < 0.05, abs(log2FoldChange) >= 0, abs(log2FoldChange) <= 1) %>% arrange(padj)
    degs_stratified$padj_buffer_05_to_10 <- res_df %>%
        filter(padj >= 0.05, padj < 0.1) %>% arrange(padj)
    degs_stratified$all_significant_padj05 <- res_df %>%
        filter(padj < 0.05) %>% arrange(padj)

    message("Stratification complete.")
    return(list(
        stratified_results = degs_stratified,
        full_results = res_df
    ))
}


# --- 3. EXAMPLE USAGE: GET "GROUND TRUTH" DEGS ---
#
message("\n--- Running 'Ground Truth' DE Analysis ---")

# --- *** NEW: SET FAST MODE HERE *** ---
# Set to a number (e.g., 50000) for a fast test run
# Set to NULL for the full production run
FAST_MODE_CELLS <- 50000
# FAST_MODE_CELLS <- NULL
# ----------------------------------------

# --- Define Parameters ---
CELL_TYPE_DIR <- "/mnt/data/shared/poweranalysis/preprocessed_data_CONSOLIDATED/Astrocyte"
DATASET_FILE <- "Roussos_Combined.qs"
SAMPLE_COL <- "donor_id"
ASSAY_NAME <- "X"
FORMULA <- ~ sex

# --- *** NEW: DYNAMIC CHECKPOINT FILE PATHS *** ---
# This adds "_temp" to filenames if in fast mode
file_suffix <- if (is.null(FAST_MODE_CELLS)) "" else "_temp"

PB_CHECKPOINT_FILE <- paste0("ground_truth_astro.pb", file_suffix, ".qs")
DISP_CHECKPOINT_FILE <- paste0("ground_truth_astro.dds.dispersions", file_suffix, ".qs")
WALD_CHECKPOINT_FILE <- paste0("ground_truth_astro.dds.wald", file_suffix, ".qs")

message(paste("Checkpoint files will be saved with suffix:", file_suffix))
# --------------------------------------------------------

# --- Load Ground Truth Data ---
if (!file.exists(PB_CHECKPOINT_FILE)) {
    message(paste("Loading full SCE file:", file.path(CELL_TYPE_DIR, DATASET_FILE)))
    ground_truth_sce <- qs::qread(file.path(CELL_TYPE_DIR, DATASET_FILE))
} else {
    message("Pseudobulk checkpoint found. Skipping initial SCE file load.")
    ground_truth_sce <- NULL
}

# --- Call the function ---
ground_truth_degs <- run_de_analysis(
    sce = ground_truth_sce,
    sample_id_col = SAMPLE_COL,
    assay_name = ASSAY_NAME,
    design_formula = FORMULA,
    pb_checkpoint = PB_CHECKPOINT_FILE,
    dds_disp_checkpoint = DISP_CHECKPOINT_FILE,
    dds_wald_checkpoint = WALD_CHECKPOINT_FILE,
    fast_subset_cells = FAST_MODE_CELLS  # <-- Pass the new parameter
)

# --- Print Summary of Ground Truth ---
if (!is.null(ground_truth_degs)) {
    message("\n--- Ground Truth Results (Adjusted P-Value Method) ---")
    message(paste("FDR < 0.05 & |logFC| > 2:",   nrow(ground_truth_degs$stratified_results$padj05_fc_gt_2)))
    message(paste("FDR < 0.05 & |logFC| 1-2:",  nrow(ground_truth_degs$stratified_results$padj05_fc_1_to_2)))
    message(paste("FDR < 0.05 & |logFC| 0-1:",  nrow(ground_truth_degs$stratified_results$padj05_fc_0_to_1)))
    message(paste("Buffer (0.05 < FDR < 0.1):",  nrow(ground_truth_degs$stratified_results$padj_buffer_05_to_10)))
    message(paste("Total Significant (FDR < 0.05):", nrow(ground_truth_degs$stratified_results$all_significant_padj05)))
}

message("\n--- Script Complete ---")
