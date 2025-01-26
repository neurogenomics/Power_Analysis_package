#' Collate DEGs detected in DGE analysis outputs, across all celltypes in a dataset (datasets/DGE analysis outputs should have common celltypes as specified below)

#' @param inpath path storing the down-sampled DGE analysis for each single-cell dataset, generated for bulk analysis
#' @param range_downsampled vector or list containing values which the data will be downsampled at, in ascending order
#' @param celltype_names list of the names specifying cell types in each DGE analysis output (in order they appear in the directory)
#' @param datasets2 list of datasets to be used, with celltype names (in DGE analysis outputs) as celltype_names2
#' @param celltype_names2 alt. list of the names specifying cell types in each DGE analysis output (in order they appear in the directory), if all are not the same as celltype
#' @param Nperms number of permutations of DGE analysis outputs for each sample
#' @param pvalue the cut-off p-value used to select DEGs

#' Saves combined list of DEGs (across all cell types) in a subdirectory inside the dataset directory

gather_celltype_DEGs <- function(inpath,
                                 range_downsampled,
                                 celltype_names,
                                 datasets2,
                                 celltype_names2,
                                 Nperms=20,
                                 pvalue=0.05){

    # validate function input params
    validate_input_parameters_bulk(inpath=inpath, range_downsampled=range_downsampled, celltype_names=celltype_names,
                                   datasets2=datasets2, celltype_names2=celltype_names2, Nperms=Nperms,
                                   pvalue=pvalue)

    # loop through datasets
    for(dataset in list.dirs(inpath,recursive=F,full.names=F)){
        print(paste0("Gathering data for ",dataset))
        # base_dataset path
        base_dataset <- paste0(inpath,"/",dataset,"/")
        # create path
        dir.create(file.path(base_dataset,"/Overall/DE_downsampling/"),recursive=TRUE,showWarnings=FALSE)
        target_base <- paste0(base_dataset,"/Overall/DE_downsampling/")
        # celltype names
        if(!dataset%in%datasets2){
            celltypes <- celltype_names
        }else{
            celltypes <- celltype_names2
        }
        # loop through sampling points
        for(sample in range_downsampled){
            # sample_samples
            sample_samples <- paste0(sample,"samples")
            print(paste0(sample_samples,"..."))
            if(dir.exists(paste0(base_dataset,"/",celltype_names[1],"/DE_downsampling/",sample_samples))){
                # loop through perms
                for(perm in 1:Nperms){
                    # list to store DEGs
                    DEGs <- c()
                    sample_perm <- paste0(sample,"_",perm)
                    # create directory in overall folder
                    dir.create(file.path(target_base,paste0(sample_samples,"/",sample_perm)),recursive=TRUE,showWarnings=FALSE)
                    target_sample_perm <- paste0(target_base,"/",sample_samples,"/",sample_perm)
                    # loop through each celltype/sample/perm
                    j <- 1
                    for(celltype in list.dirs(base_dataset,recursive=F,full.names=F)){
                        if(celltype!="Overall"){
                            # go inside celltype/sample/perm
                            cell_sample_permpath <- paste0(base_dataset,"/",celltype,"/DE_downsampling/",sample_samples,"/",sample_perm)
                            # load DGE analysis output
                            load(paste0(cell_sample_permpath,"/DEout",sample_perm,".RData"))
                            # get DEGs
                            celltype_DE <- celltypes[[j]]
                            DEGs <- c(DEGs, subset(eval(as.name(paste0("DEout_",sample)))$celltype_all_genes[[celltype_DE]],PValue<=pvalue)$name)
                            j <- j+1
                        }
                    }
                    DEGs <- unique(DEGs)
                    # save DEGs
                    save(DEGs, file=paste0(target_sample_perm,"/DEout",sample_perm,".RData"))
                }
            }
        }
    }
}
