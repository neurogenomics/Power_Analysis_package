# SCRIPT 1: DATA PRE-PROCESSING (Re-structured Output)
#
# GOAL:
# 1. Load each raw dataset one-by-one (sequentially).
# 2. Identify the cell type column.
# 3. For each cell type, save the subsetted data to a top-level
#    cell type folder, named after the dataset.
# 4. Clear the large dataset from memory before loading the next.
#
# NEW OUTPUT STRUCTURE:
#   /preprocessed_data/
#   |-- Astro/
#   |   |-- Roussos_Combined.qs
#   |   |-- Allen_Combined.qs
#   |   |-- Tsai.qs
#   |-- Microglia/
#   |   |-- Roussos_Combined.qs
#   |   |-- Allen_Combined.qs
#   ...

# --- 1. SETUP ---
# Load necessary libraries
library(SingleCellExperiment)
library(zellkonverter)
library(qs)

# --- Define Paths ---
# INPUT: Path to the raw .h5ad and .qs files
shared.folder <- "/mnt/data/shared/poweranalysis/raw_data"

# OUTPUT: Path to the new pre-processed data store
output.base.dir <- "/mnt/data/shared/poweranalysis/preprocessed_data"

# Create the main output directory if it doesn't exist
if (!dir.exists(output.base.dir)) {
    dir.create(output.base.dir, recursive = TRUE)
}

# --- 2. HELPER FUNCTION (MODIFIED) ---

#' Processes a full SCE object and saves it in a structured directory.
#'
#' @param sce A loaded SingleCellExperiment object.
#' @param dataset.name The name to use for the main dataset folder (e.g., "Roussos_Combined").
#' @param celltype.col The string name of the column in colData(sce) that
#'                     contains the cell type labels (e.g., "class" or "cell_type").
#' @param output.base.dir The root directory to save the data to.
process_and_save_dataset <- function(sce, dataset.name, celltype.col, output.base.dir) {
    message(paste0("\n--- Processing: ", dataset.name, " ---"))

    # Get the list of unique cell types, removing any NAs
    cell.types <- try(unique(na.omit(sce[[celltype.col]])), silent = TRUE)

    # Handle cases where the celltype.col doesn't exist
    if (inherits(cell.types, "try-error")) {
        warning(paste("Could not find cell type column:", celltype.col, "in dataset:", dataset.name))
        return()
    }

    cell.types <- cell.types[cell.types != ""] # Remove empty strings if any

    if (length(cell.types) == 0) {
        warning(paste("No cell types found for", dataset.name, "using column", celltype.col))
        return()
    }

    message(paste("Found", length(cell.types), "cell types:", paste(cell.types, collapse=", ")))

    # Loop through each cell type, subset, and save
    for (cell in cell.types) {

        # Sanitize cell type name for directory path
        cell.name.safe <- gsub("[^a-zA-Z0-9_.-]", "_", cell)

        # --- NEW DIRECTORY LOGIC ---
        # 1. Create the top-level cell type folder
        cell.dir <- file.path(output.base.dir, cell.name.safe)
        if (!dir.exists(cell.dir)) {
            dir.create(cell.dir, recursive = TRUE)
            message(paste("Created new cell type directory:", cell.dir))
        }

        # 2. Define the output file path using the dataset name
        output.file <- file.path(cell.dir, paste0(dataset.name, ".qs"))
        # --- END NEW LOGIC ---

        # Subset the dataset
        message(paste("Subsetting for:", cell))
        # Note: Using sce[[celltype.col]] syntax for subsetting
        subset.sce <- subset(sce, , sce[[celltype.col]] == cell)

        # Save the subsetted object
        message(paste("Saving", dataset.name, "data to:", output.file))
        qs::qsave(subset.sce, output.file)

        # Clean up subset
        rm(subset.sce)
    }

    message(paste0("--- Finished processing: ", dataset.name, " ---"))
}


# --- 3. SEQUENTIAL EXECUTION ---
#
# We will load one "group" of data, process it, save it,
# then remove it from memory (rm() and gc()) before starting the next.

# ---
# Dataset 1: Roussos Combined (PsychAD)
# ---
try({
    print("Loading Roussos 1 (MSSM)...")
    roussos1 <- zellkonverter::readH5AD(file.path(shared.folder, "PsychAD-MSSM-e6600c0c-5930-4ca6-b872-23a4fb33d9e4.h5ad"), use_hdf5 = TRUE)
    print("Loading Roussos 2 (HBCC)...")
    roussos2 <- zellkonverter::readH5AD(file.path(shared.folder, "PsychAD-HBCC-00056957-c963-493a-aaff-2a5fa3047a06.h5ad"), use_hdf5 = TRUE)
    print("Loading Roussos 3 (RADC)...")
    roussos3 <- zellkonverter::readH5AD(file.path(shared.folder, "PsychAD-RADC-a4bc87dd-59a2-4bdf-bb88-da2c1ff6daf5.h5ad"), use_hdf5 = TRUE)

    print("Combining Roussos datasets...")
    roussos.combined <- SingleCellExperiment::cbind(roussos1, roussos2, roussos3)

    print("Cleaning up individual Roussos objects...")
    rm(roussos1, roussos2, roussos3)
    gc() # Garbage collect

    # Process the combined object. From the original script, the cell type column is "class"
    process_and_save_dataset(roussos.combined, "Roussos_Combined", "class", output.base.dir)

    print("Cleaning up combined Roussos object...")
    rm(roussos.combined)
    gc() # Garbage collect

}, silent = FALSE)


# ---
# Dataset 2: Allen Combined
# ---
try({
    print("Loading Allen datasets...")
    allen_astro <- qs::qread(file.path(shared.folder, "allen_astrocytes.qs"))
    allen_endo <- qs::qread(file.path(shared.folder, "allen_Endo.qs"))
    allen_oligo <- qs::qread(file.path(shared.folder, "allen_oligo.qs"))
    allen_micro <- qs::qread(file.path(shared.folder, "allen_microglia.qs"))

    # Pre-processing from original script
    colnames(colData(allen_astro)) <- sub("Donor.ID", "donor_id", colnames(colData(allen_astro)))
    colnames(colData(allen_micro)) <- sub("Donor.ID", "donor_id", colnames(colData(allen_micro)))

    common_cols <- intersect(intersect(
        intersect(colnames(colData(allen_astro)), colnames(colData(allen_endo))),
        colnames(colData(allen_micro))),
        colnames(colData(allen_oligo)))

    colData(allen_astro) <- colData(allen_astro)[, common_cols]
    colData(allen_endo) <- colData(allen_endo)[, common_cols]
    colData(allen_micro) <- colData(allen_micro)[, common_cols]
    colData(allen_oligo) <- colData(allen_oligo)[, common_cols]

    print("Combining Allen datasets...")
    allen.combined <- SingleCellExperiment::cbind(allen_astro, allen_endo, allen_micro, allen_oligo)

    print("Cleaning up individual Allen objects...")
    rm(allen_astro, allen_endo, allen_oligo, allen_micro)
    gc()

    # Process. From original script, cell type column is "cell_type"
    process_and_save_dataset(allen.combined, "Allen_Combined", "cell_type", output.base.dir)

    print("Cleaning up combined Allen object...")
    rm(allen.combined)
    gc()

}, silent = FALSE)


# ---
# Other Datasets (Processed one-by-one)
#
# Based on the `celltypeIDs_ind` vector in the original script:
# "class" (Roussos)
# "cell_type" (Allen)
# "cluster_celltype" (Gerrits_EC)
# (missing) (Gerrits_OC) -> Assume "cluster_celltype"
# "cluster_celltype" (Gerrits_OTC)
# "cluster_celltype" (Gerrits_SSC)
# "cluster_celltype" (ICL_MTG)
"cluster_celltype" (ICL_SSC)
# "cluster_celltype" (Smith_EC)
# "cluster_celltype" (Smith_SSC)
# "cluster_celltype" (Tsai)
# "cluster_celltype" (Zhou)
#
# It seems "cluster_celltype" is the default for most.
# ---

# List of remaining datasets to process
# Format: list(list(file_name, dataset_name, celltype_col), ...)
datasets_to_process <- list(
    list("Tsai.qs", "Tsai", "cluster_celltype"),
    list("Zhou.qs", "Zhou", "cluster_celltype"),
    list("Gerrits_EC.qs", "Gerrits_EC", "cluster_celltype"),
    list("Gerrits_OC.qs", "Gerrits_OC", "cluster_celltype"),
    list("Gerrits_OTC.qs", "Gerrits_OTC", "cluster_celltype"),
    list("Gerrits_SSC.qs", "Gerrits_SSC", "cluster_celltype"),
    list("Trem2_MTG.qs", "ICL_Cortex_MTG", "cluster_celltype"),
    list("Trem2_SSC.qs", "ICL_Cortex_SSC", "cluster_celltype"),
    list("Smith_EC.qs", "Smith_EC", "cluster_celltype"),
    list("Smith_SSC.qs", "Smith_SSC", "cluster_celltype")
)

# Loop and process each remaining dataset
for (ds in datasets_to_process) {
    file.name <- ds[[1]]
    dataset.name <- ds[[2]]
    celltype.col <- ds[[3]]

    try({
        print(paste("Loading dataset:", dataset.name))
        sce.obj <- qs::qread(file.path(shared.folder, file.name))

        process_and_save_dataset(sce.obj, dataset.name, celltype.col, output.base.dir)

        print(paste("Cleaning up dataset:", dataset.name))
        rm(sce.obj)
        gc()

    }, silent = FALSE)
}

message("\n--- ALL DATA PRE-PROCESSING COMPLETE ---")
message(paste("All subsetted data has been saved to:", output.base.dir))
