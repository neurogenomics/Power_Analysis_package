# Define global variables
utils::globalVariables(c("PValue","name"))

#' Downsample the dataset, based either on the individuals or cells, and run DE analysis on each downsampled output. Save results in a dataframe

#' @importFrom stats as.formula

#' @param data the input data (should be an SCE object)
#' @param range_downsampled range of values to be downsampled for, in ascending order
#' @param output_path base path in which outputs will be stored
#' @param sampled downsampling carried out based on what (either "individuals" or "cells")
#' @param sampleID sample ID
#' @param design the design formula of class type `formula`. Equation used to fit the model- data for the generalised linear model e.g. expression ~ sex + pmi + disease
#' @param sexID sex ID
#' @param celltypeID cell type ID
#' @param coeff which coefficient to carry out DE analysis with respect to
#' @param fdr the cut-off False Discovery Rate below which to select DEGs
#' @param nom_pval the cut-off nominal P-value below which to select DEGs (as an alternative to FDR)
#' @param Nperms number of subsets created when downsampling at each level
#' @param y the column name in the SCE object for the return variable e.g. "diagnosis" - Case or disease. Default is the last variable in the design formula. y can be discrete (logistic regression) or continuous (linear regression)
#' @param region the column name in the SCE object for the study region. If there are multiple regions in the study (for example two brain regions). Pseudobulk values can be derived separately. Default is "single_region" which will not split by region.
#' @param control character specifying which control level for the differential expression analysis e.g. in a case/control/other study use "control" in the y column to compare against. NOTE only need to specify if more than two groups in y, leave as default value for two groups or continuous y. Default is NULL.
#' @param pval_adjust_method the adjustment method for the p-value in the differential expression analysis. Default is benjamini hochberg "BH". See  stats::p.adjust for available options
#' @param rmv_zero_count_genes whether genes with no count values in any cell should be removed. Default is TRUE

#' Saves all DE outputs for downsampled files as well as a summary table of results showing number of true DEGs detected at each number of samples/cells

downsampling_DEanalysis <- function(data,
                                    range_downsampled="placeholder",
                                    output_path=getwd(),
                                    sampled="individuals",
                                    sampleID="donor_id",
                                    design="placeholder",
                                    sexID="sex",
                                    celltypeID="cell_type",
                                    coeff="male",
                                    fdr=0.05,
                                    nom_pval=0.05,
                                    Nperms=20,
                                    y=NULL,
                                    region="single_region",
                                    control=NULL,
                                    pval_adjust_method="BH",
                                    rmv_zero_count_genes=TRUE){

    # alter range_downsampled
    if(range_downsampled=="placeholder"){
        range_downsampled <- downsampling_range(data, sampled, sampleID)
    }
    # alter design
    if(design=="placeholder"){
        design=as.formula(paste0("~",sexID))
    }

    # create output path if doesn't already exist
    dir.create(output_path,showWarnings=FALSE)
    setwd(output_path)
    # validate function input params
    validate_input_parameters_power(data=data, range_downsampled=range_downsampled, output_path=output_path,
                                    sampled=sampled, sampleID=sampleID, design=design, 
                                    sexID=sexID, celltypeID=celltypeID, coeff=coeff,
                                    fdr=fdr, nom_pval=nom_pval, Nperms=Nperms,
                                    y=y, region=region, control=control,
                                    pval_adjust_method=pval_adjust_method, rmv_zero_count_genes=rmv_zero_count_genes)

    # get celltype name from dataset
    celltype_name <- toString(unique(data[[celltypeID]]))
    # check if DE analysis output present already in output_path
    if(!"DEout.RData" %in% list.files(output_path)){
        # run and save DE analysis
        assign("DEout", sc_cell_type_de(data, design=design, pseudobulk_ID=sampleID, celltype_ID=celltypeID, y=y, region=region, control=control, pval_adjust_method=pval_adjust_method, rmv_zero_count_genes=rmv_zero_count_genes, verbose=T, coef=coeff))
        save(DEout,file=paste0(output_path,"/DEout.RData"))
    }else{
        load(paste0(output_path,"/DEout.RData"))
    }
    # get DEGs using both FDR and nominal Pval
    DEGs_fdr <- subset(DEout$celltype_all_genes[[celltype_name]], adj_pval<fdr)$name
    DEGs_pval <- subset(DEout$celltype_all_genes[[celltype_name]], PValue<nom_pval)$name
	# all genes
	allgenes <- DEout$celltype_all_genes[[celltype_name]]
	# non DEGs
	nonDEGs_fdr <- subset(allgenes,!name%in%DEGs_fdr)$name
	nonDEGs_pval <- subset(allgenes,!name%in%DEGs_pval)$name
    # create dfs
    Iteration <- 1:Nperms
    DEGs_detected_fdr <- data.frame(Iteration)
    DEGs_detected_pval <- data.frame(Iteration)
	falsePositives_fdr <- data.frame(Iteration)
	falsePositives_pval <- data.frame(Iteration)

    ### downsampling, DGE analysis
    if(sampled=="individuals"){

        # make sure directory is correct (to save output)
        output_path <- paste0(output_path,"/DE_downsampling/")
        # create path
        dir.create(output_path,showWarnings=FALSE)
        # loop through downsampled values
        for(value in range_downsampled){
            # create and redirect to new folder
            path_val <- paste0(output_path,"/",toString(value),"samples")
            dir.create(path_val,showWarnings=FALSE)
            setwd(path_val)
            # create subsets
            samples <- sample_individuals(data, value, sampleID, sexID, Nperms)
            # run DGE analysis on each subset, get number of DEGs and FPs
            assign(paste0("numDEGs_fdr_",toString(value)), list())
            assign(paste0("numDEGs_pval_",toString(value)), list())
			assign(paste0("numFPs_fdr_",toString(value)), list())
            assign(paste0("numFPs_pval_",toString(value)), list())
            # loop through
            for(j in 1:Nperms){
                # create sub-directory for each one
                path <- paste0(path_val,"/",toString(value),"_",j) 
                dir.create(path,showWarnings=FALSE)
                setwd(path)
                # ensure sexID isnt "Sex", has to be lower case (change this in SCE if needed)
                assign(paste0("DEout_",toString(value)), sc_cell_type_de(samples[[j]], design=design, pseudobulk_ID=sampleID, celltype_ID=celltypeID, y=y, region=region, control=control, pval_adjust_method=pval_adjust_method, rmv_zero_count_genes=rmv_zero_count_genes, verbose=T, coef=coeff))
                # save output
                save(list=eval(paste0("DEout_",toString(value))),file=paste0("DEout",toString(value),"_",j,".RData"))
                # get number of TP DEGs
                degs_fdr <- subset(eval(as.name(paste0("DEout_",toString(value))))$celltype_all_genes[[celltype_name]], adj_pval<fdr)$name
                degs_pval <- subset(eval(as.name(paste0("DEout_",toString(value))))$celltype_all_genes[[celltype_name]], PValue<nom_pval)$name
                eval(parse(text=paste0(paste0("numDEGs_fdr_",toString(value)), "[[j]]<-", length(degs_fdr[degs_fdr %in% DEGs_fdr]))))
                eval(parse(text=paste0(paste0("numDEGs_pval_",toString(value)), "[[j]]<-", length(degs_pval[degs_pval %in% DEGs_pval]))))
				# get number of FP DEGs
				eval(parse(text=paste0(paste0("numFPs_fdr_",toString(value)), "[[j]]<-", length(degs_fdr[degs_fdr %in% nonDEGs_fdr]))))
                eval(parse(text=paste0(paste0("numFPs_pval_",toString(value)), "[[j]]<-", length(degs_pval[degs_pval %in% nonDEGs_pval]))))
            }
            # unlist            
            assign(paste0("numDEGs_fdr_",toString(value)), unlist(eval(as.name(paste0("numDEGs_fdr_",toString(value))))))
            assign(paste0("numDEGs_pval_",toString(value)), unlist(eval(as.name(paste0("numDEGs_pval_",toString(value))))))
			assign(paste0("numFPs_fdr_",toString(value)), unlist(eval(as.name(paste0("numFPs_fdr_",toString(value))))))
            assign(paste0("numFPs_pval_",toString(value)), unlist(eval(as.name(paste0("numFPs_pval_",toString(value))))))
            # add to dataframe
            DEGs_detected_fdr[[paste0("numDEGs_",toString(value),"samples")]] <- eval(as.name(paste0("numDEGs_fdr_",toString(value))))
            DEGs_detected_pval[[paste0("numDEGs_",toString(value),"samples")]] <- eval(as.name(paste0("numDEGs_pval_",toString(value))))
			falsePositives_fdr[[paste0("numFPs_",toString(value),"samples")]] <- eval(as.name(paste0("numFPs_fdr_",toString(value))))
			falsePositives_pval[[paste0("numFPs_",toString(value),"samples")]] <- eval(as.name(paste0("numFPs_pval_",toString(value))))
        }
        setwd(output_path)
        save(DEGs_detected_fdr, file="DEGs_detected_fdr.RData")
        save(DEGs_detected_pval, file="DEGs_detected_pval.RData")
		save(falsePositives_fdr, file="falsePositives_fdr.RData")
		save(falsePositives_pval, file="falsePositives_pval.RData")

    }else{

        # make sure directory is correct (to save output)
        output_path <- paste0(output_path,"/DE_downsampling_cells/")
        # create path
        dir.create(output_path,showWarnings=FALSE)
        # loop through downsampled values
        for(value in range_downsampled){
            # create and redirect to new folder
            path_val <- paste0(output_path,"/",toString(value),"cells_persample")
            dir.create(path_val,showWarnings=FALSE)
            setwd(path_val)
            # create subsets
            cells <- sample_cells(data, value, sampleID, Nperms)
            # run DGE analysis on each subset, get number of DEGs
            assign(paste0("numDEGs_fdr_",toString(value)), list())
            assign(paste0("numDEGs_pval_",toString(value)), list())
			assign(paste0("numFPs_fdr_",toString(value)), list())
            assign(paste0("numFPs_pval_",toString(value)), list())			
            # loop through
            for(j in 1:Nperms){
                # create sub-directory for each one
                path <- paste0(path_val,"/",toString(value),"_",j) 
                dir.create(path,showWarnings=FALSE)
                setwd(path)
                # ensure sexID isnt "Sex", has to be lower case (change this in SCE if needed)
                assign(paste0("DEout_",toString(value)), sc_cell_type_de(cells[[j]], design=design, pseudobulk_ID=sampleID, celltype_ID=celltypeID, y=y, region=region, control=control, pval_adjust_method=pval_adjust_method, rmv_zero_count_genes=rmv_zero_count_genes, verbose=T, coef=coeff))
                # save output
                save(list=eval(paste0("DEout_",toString(value))),file=paste0("DEout",toString(value),"_",j,".RData"))
                # get number of TP DEGs
                degs_fdr <- subset(eval(as.name(paste0("DEout_",toString(value))))$celltype_all_genes[[celltype_name]], adj_pval<fdr)$name
                degs_pval <- subset(eval(as.name(paste0("DEout_",toString(value))))$celltype_all_genes[[celltype_name]], PValue<nom_pval)$name
                eval(parse(text=paste0(paste0("numDEGs_fdr_",toString(value)), "[[j]]<-", length(degs_fdr[degs_fdr %in% DEGs_fdr]))))
                eval(parse(text=paste0(paste0("numDEGs_pval_",toString(value)), "[[j]]<-", length(degs_pval[degs_pval %in% DEGs_pval]))))
				# get number of FP DEGs
				eval(parse(text=paste0(paste0("numFPs_fdr_",toString(value)), "[[j]]<-", length(degs_fdr[degs_fdr %in% nonDEGs_fdr]))))
                eval(parse(text=paste0(paste0("numFPs_pval_",toString(value)), "[[j]]<-", length(degs_pval[degs_pval %in% nonDEGs_pval]))))
            }
            # unlist            
            assign(paste0("numDEGs_fdr_",toString(value)), unlist(eval(as.name(paste0("numDEGs_fdr_",toString(value))))))
            assign(paste0("numDEGs_pval_",toString(value)), unlist(eval(as.name(paste0("numDEGs_pval_",toString(value))))))
			assign(paste0("numFPs_fdr_",toString(value)), unlist(eval(as.name(paste0("numFPs_fdr_",toString(value))))))
            assign(paste0("numFPs_pval_",toString(value)), unlist(eval(as.name(paste0("numFPs_pval_",toString(value))))))			
            # add to dataframe
            DEGs_detected_fdr[[paste0("numDEGs_",toString(value),"cells_persample")]] <- eval(as.name(paste0("numDEGs_fdr_",toString(value))))
            DEGs_detected_pval[[paste0("numDEGs_",toString(value),"cells_persample")]] <- eval(as.name(paste0("numDEGs_pval_",toString(value))))
			falsePositives_fdr[[paste0("numFPs_",toString(value),"samples")]] <- eval(as.name(paste0("numFPs_fdr_",toString(value))))
			falsePositives_pval[[paste0("numFPs_",toString(value),"samples")]] <- eval(as.name(paste0("numFPs_pval_",toString(value))))			
        }
        setwd(output_path)
        save(DEGs_detected_fdr, file="DEGs_detected_fdr.RData")
        save(DEGs_detected_pval, file="DEGs_detected_pval.RData")
		save(falsePositives_fdr, file="falsePositives_fdr.RData")
		save(falsePositives_pval, file="falsePositives_pval.RData")		

    }
}