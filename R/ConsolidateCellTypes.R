# SCRIPT 2: CONSOLIDATE AND HARMONIZE CELL TYPES
#
# GOAL:
# 1. Take the output from script 1 (which is fragmented by
#    dataset-specific cell names).
# 2. Define a "consolidation map" to group similar cell types.
# 3. Create a NEW, clean directory structure where all files for
#    a harmonized cell type (e.g., "Astrocyte") are in one folder.

# --- 1. SETUP ---

# INPUT: The output directory from script 1
# (Change this if it's different)
input.dir <- "/mnt/data/shared/poweranalysis/preprocessed_data"

# OUTPUT: The new, clean, consolidated directory
output.dir <- "/mnt/data/shared/poweranalysis/preprocessed_data_CONSOLIDATED"

if (!dir.exists(output.dir)) {
    dir.create(output.dir, recursive = TRUE)
}

# --- 2. DEFINE THE CONSOLIDATION MAP ---
# This list defines our "merge" rules.
#
# Format: "New_Folder_Name" = c("old_folder_name_1", "old_folder_name_2", ...)
#
# We are mapping all the fragmented names from your list to a
# single, harmonized name.

consolidation_map <- list(

    # Astrocytes
    "Astrocyte" = c(
        "Astro",
        "astrocyte_of_the_cerebral_cortex"
    ),

    # Excitatory Neurons
    "Excitatory_Neuron" = c(
        "EN",
        "EN-L2", "EN-L2-3", "EN-L2-4",
        "EN-L3-4", "EN-L3-5",
        "EN-L4-5", "EN-L4-6",
        "EN-L5-6", "EN-L5-6-I", "EN-L5-6-II"
    ),

    # Inhibitory Neurons
    "Inhibitory_Neuron" = c(
        "IN",
        "IN-GAD1",
        "IN-LAMP5",
        "IN-PVALB",
        "IN-SST", "IN-SST-NMBR", "IN-SST-NPY",
        "IN-VIP"
    ),

    # Endothelial
    "Endothelial" = c(
        "Endo",
        "cerebral_cortex_endothelial_cell"
    ),

    # Microglia / Immune
    "Microglia_Immune" = c(
        "Immune",
        "Micro",
        "microglial_cell"
    ),

    # Oligodendrocytes
    "Oligodendrocyte" = c(
        "Oligo",
        "oligodendrocyte"
    ),

    # OPCs
    "OPC" = c(
        "OPC"
    ),

    # Mural Cells
    "Mural" = c(
        "Mural",
        "Pericyte",
        "VSMC"
    )
)

# --- 3. RUN THE CONSOLIDATION ---

message("Starting consolidation...")
message(paste("Source:", input.dir))
message(paste("Destination:", output.dir))

# Loop over each NEW, clean folder name
for (new_folder_name in names(consolidation_map)) {

    message(paste0("\n--- Processing: ", new_folder_name, " ---"))

    # 1. Create the new destination folder
    new_dest_dir <- file.path(output.dir, new_folder_name)
    if (!dir.exists(new_dest_dir)) {
        dir.create(new_dest_dir, recursive = TRUE)
    }

    # 2. Get the list of old, fragmented folders to merge
    source_folders <- consolidation_map[[new_folder_name]]

    # 3. Loop over the old folders and copy their contents
    for (old_folder_name in source_folders) {

        old_source_dir <- file.path(input.dir, old_folder_name)

        # Check if the old source directory actually exists
        if (dir.exists(old_source_dir)) {

            # List all the .qs files in it
            files_to_copy <- list.files(
                old_source_dir,
                pattern = "\\.qs$",
                full.names = TRUE
            )

            if (length(files_to_copy) > 0) {
                message(paste("Copying files from:", old_source_dir))

                # Copy them to the new consolidated folder
                file.copy(from = files_to_copy, to = new_dest_dir)

            } else {
                message(paste("No .qs files found in:", old_source_dir))
            }

        } else {
            warning(paste("Source folder not found, skipping:", old_source_dir))
        }
    }
}

message("\n--- CONSOLIDATION COMPLETE ---")
message(paste("All files have been merged into:", output.dir))
