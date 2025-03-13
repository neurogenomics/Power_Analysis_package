#' Collate DEGs detected in DGE analysis outputs, across all celltypes in a dataset (datasets/DGE analysis outputs should have common celltypes as specified below)

#' @param range_downsampled vector or list containing values which the data will be downsampled at, in ascending order
#' @param celltype_correspondence list of different names specifying each cell type
#' @param Nperms number of permutations of DGE analysis outputs for each sample
#' @param pvalue the cut-off p-value used to select DEGs
#' @param output_path path storing the down-sampled DGE analysis for each single-cell dataset, generated for bulk analysis

#' Saves combined list of DEGs (across all cell types) in a subdirectory inside the dataset directory

gather_celltype_DEGs <- function(range_downsampled,
                                 celltype_correspondence,
                                 Nperms=20,
                                 pvalue=0.05,
                                 output_path=getwd()){

    # validate function input params
    validate_input_parameters_bulk(output_path=output_path, range_downsampled=range_downsampled, celltype_correspondence=celltype_correspondence,
                                   Nperms=Nperms, pvalue=pvalue)

    # loop through datasets
    for(dataset in list.dirs(output_path,recursive=F,full.names=F)){
        print(paste0("Gathering data for ",dataset))
        # base_dataset path
        base_dataset <- file.path(output_path,dataset)
        # create path
        target_base <- file.path(base_dataset,"Overall/DE_downsampling")
        dir.create(target_base,recursive=TRUE,showWarnings=FALSE)
        # loop through sampling points
        for(sample in range_downsampled){
            # sample_samples
            sample_samples <- paste0(sample,"samples")
            print(paste0(sample_samples,"..."))
            # loop through perms
            for(perm in 1:Nperms){
                # list to store DEGs
                DEGs <- c()
                sample_perm <- paste0(sample,"_",perm)
                # create directory in overall folder
                dir.create(file.path(target_base,sample_samples,sample_perm),recursive=TRUE,showWarnings=FALSE)
                target_sample_perm <- file.path(target_base,sample_samples,sample_perm)
                # loop through each celltype
                for(standard_celltype in names(celltype_correspondence)){
                    celltype_names <- celltype_correspondence[[standard_celltype]]
                    # check if any of the cell type directories exist
                    celltype_dir <- NULL
                    for (celltype_name in celltype_names) {
                        if (dir.exists(file.path(base_dataset, celltype_name, "DE_downsampling"))) {
                        celltype_dir <- celltype_name
                        break
                        }
                    }
                    if (is.null(celltype_dir)) {
                        next
                    }
                    cell_sample_permpath <- file.path(base_dataset,celltype_dir,"DE_downsampling",sample_samples,sample_perm)
                    # load DGE analysis output
                    deg_file <- file.path(cell_sample_permpath,"DEout",paste0(sample_perm,".RData"))
                    if(file.exists(deg_file)){
                        load(deg_file)
                        # get DEGs
                        DEGs <- c(DEGs, subset(eval(as.name(paste0("DEout_",sample)))$celltype_all_genes[[celltype_dir]],PValue<=pvalue)$name)
                        j <- j+1
                        DEGs <- unique(DEGs)
                    }
                }
                # save DEGs
                save(DEGs, file=file.path(target_sample_perm,"DEout",paste0(sample_perm,".RData")))
            }
        }
    }
}
