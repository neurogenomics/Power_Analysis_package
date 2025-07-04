# Define global variables
utils::globalVariables(c("deg_direction",".I","adj_pval","celltype","i.deg_direction","phenotype","gene_name","x.hgnc_symbol",".N","N","num_cells","prop","N_prop","colour_ident","."))

#' Create differential expression analysis plots. Run by DGE_analysis()

#' @importFrom EnsDb.Hsapiens.v79 EnsDb.Hsapiens.v79
#' @importFrom data.table rbindlist setnames as.data.table setkey data.table setorder :=
#' @importFrom ggplot2 ggplot geom_jitter stat_summary scale_shape_manual labs facet_wrap theme ggsave geom_bar ggtitle geom_hline scale_colour_manual aes element_text geom_point element_blank unit geom_text mean_cl_normal
#' @importFrom reshape2 melt
#' @importFrom ensembldb select
#' @importFrom cowplot theme_cowplot
#' @importFrom viridis scale_colour_viridis scale_fill_viridis
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats median quantile
#' @importFrom utils globalVariables
#' @importFrom scales percent
#' @importFrom grid unit

#' @param pb_dat A list containing
#'     sumDat: matrix of the summed pseudobulk count values
#'     annot_pb: dataframe of the annotation data from the SCE rolled up based on the pseudobulk aggregation.
#' @param y the column name in the SCE object for the return variable e.g. "diagnosis" - Case or disease. y can be discrete (logisitc regression) or continuous (linear regression)
#' @param celltype_DEGs_dt data table containing the DEGs for each cell type with their differential expression data
#' @param celltype_all_genes_dt data table containing the all genes for each cell type with their differential expression data
#' @param counts_celltypes vector with the counts of cells after QC in each cell type
#' @param output_path the folder where the graphs from the differential expression analysis are saved. Default will create a folder in the current working directory "sc.cell.type.de.graphs". False will skip plotting.

#' Creates and saves plot in relevant directory

plot_de_analysis <- function(pb_dat,y,celltype_DEGs_dt,celltype_all_genes_dt,
                             counts_celltypes,output_path){
    # save(pb_dat,y,celltype_DEGs_dt,celltype_all_genes_dt,
    #      counts_celltypes,output_path, file="temp.RData")
    logFC = name = NULL
    #if it doesn't exist already make folder for plots
    dir.create(output_path,showWarnings = F,recursive=TRUE)

    celltype_DEGs_dt[,deg_direction:="Down"]
    celltype_DEGs_dt[logFC>0,deg_direction:="Up"]
    #get top 3 most sig p-values for each cell type
    top_up_down_genes <-celltype_DEGs_dt[
        celltype_DEGs_dt[,.I[adj_pval %in% sort(adj_pval)[1:3]],by=celltype]$V1]
    #plot pseudobulk expression
    #first filter cell type expression to genes of interest
    top_degs_pseudobulk_exp <- pb_dat[unique(top_up_down_genes$celltype)]
    for(ct_i in names(top_degs_pseudobulk_exp)){
        ct_genes <- top_up_down_genes[celltype==ct_i,name]
        #filter and melt to get long DF of genes of interest
        #/sum(pb_dat$Micro$sumDat)
        top_degs_pseudobulk_exp[[ct_i]] <-
            melt(top_degs_pseudobulk_exp[[ct_i]]$sumDat[ct_genes,,drop=FALSE])
    }
    #combine
    top_degs_pseudobulk_exp <-
        rbindlist(top_degs_pseudobulk_exp,idcol="celltype")
    setnames(top_degs_pseudobulk_exp,c("celltype","name","group_sample",
                                        "expression"))
    #add on DEG direction
    top_degs_pseudobulk_exp[top_up_down_genes,deg_direction:=i.deg_direction,
                                on=c("name","celltype")]
    #add in y - first get annotation data in single datatable
    annot_dt <- lapply(pb_dat, function(x) x$annot_pb)
    annot_dt <- rbindlist(annot_dt,idcol="celltype")
    top_degs_pseudobulk_exp[annot_dt,phenotype:=get(y),
                            on=c("group_sample","celltype")]
    #add in gene names
    genes <- unique(as.character(top_degs_pseudobulk_exp$name))
    gene_IDs <- select(EnsDb.Hsapiens.v79, keys= genes, keytype = "GENEID", columns = c("GENEID","SYMBOL"))
    colnames(gene_IDs) <- c("ensembl_gene_id","hgnc_symbol")
    gene_IDs <- as.data.table(gene_IDs)
    setnames(gene_IDs,"ensembl_gene_id","name")
    #remove any dups in the reference set - two names for one ENSEMBL ID
    gene_IDs <- unique(gene_IDs,by="name")
    setkey(top_degs_pseudobulk_exp,name)
    #append gene names
    top_degs_pseudobulk_exp[, gene_name := gene_IDs
                            [top_degs_pseudobulk_exp, on=.(name),
                                    x.hgnc_symbol]]
    #plot increase size A4
    top_degs_pseudobulk_exp_plot <-
        ggplot(top_degs_pseudobulk_exp[gene_name!=""&!is.na(deg_direction),],
               aes(x = phenotype, y = expression,colour=deg_direction)) +
        geom_jitter(height=0) +
        stat_summary(fun.data = "mean_cl_normal",
                     colour = "grey",
                     geom = "crossbar",
                     show.legend = T)+
        labs(y= "Sum pseudobulk expression counts (unnormalised)",
                x = "Phenotype",
             colour="DEG direction") +
        facet_wrap(~ celltype+gene_name, scales = "free_y")+
        theme_cowplot()+
        theme(strip.text.x = element_text(size = 6),
              axis.text = element_text(size=6))+
        scale_colour_viridis(discrete = T)
    #save the graph to output_path
    suppressMessages(ggsave(path = output_path,
                            filename = "Pseudobulk_exp_most_sig_genes.pdf",
                            plot=top_degs_pseudobulk_exp_plot,
                            width = 210, height = 297, units = "mm"))
    suppressMessages(ggsave(path = output_path,
                            filename = "Pseudobulk_exp_most_sig_genes.png",
                            plot=top_degs_pseudobulk_exp_plot,
                            width = 210, height = 297, units = "mm"))

    #plot number of cells per cell type
    counts_celltypes_dt <- data.table(celltype=names(counts_celltypes),
                                      counts=counts_celltypes)
    cell_counts_plot <-
        ggplot(data=counts_celltypes_dt,
               aes(x=factor(celltype),y = counts, fill=factor(celltype)))+
        geom_bar(stat="identity")+
        labs(y= "Number of cells after QC", x = "Cell Type",fill="Cell Type") +
        geom_text(aes(label = counts), vjust = -0.25, size = 3) +
        theme_cowplot()+
        theme(axis.text = element_text(size=6))+
        scale_fill_viridis(discrete = T)
    #save the graph to output_path
    suppressMessages(ggsave(path = output_path,
                            filename = "Cell_counts_after_QC.pdf",
                            plot=cell_counts_plot))
    suppressMessages(ggsave(path = output_path,
                            filename = "Cell_counts_after_QC.png",
                            plot=cell_counts_plot))

    #plot DEGs per cell type
    deg_per_cell_type_plot <-
        ggplot(data = celltype_DEGs_dt[,.N,by=celltype],
               aes(x = factor(celltype), y = N, fill = factor(celltype))) +
        geom_bar(stat = "identity") +
        labs(y = "Number of DEGs identified", x = "Cell Type", fill = "Cell Type") +
        geom_text(aes(label = N), vjust = -0.25, size = 3) +
        theme_cowplot() +
        theme(axis.text = element_text(size = 9)) +
        scale_fill_viridis(discrete = T)
    #save the graph to output_path
    suppressMessages(ggsave(path = output_path,
                            filename = "DEGs_per_cell_type.pdf",
                            plot=deg_per_cell_type_plot))
    suppressMessages(ggsave(path = output_path,
                            filename = "DEGs_per_cell_type.png",
                            plot=deg_per_cell_type_plot))

    #plot degs as proportion of cell numbers
    degs_prop <- celltype_DEGs_dt[,.N,by=celltype]
    setkey(degs_prop,celltype)
    degs_prop[names(counts_celltypes),num_cells:=counts_celltypes]
    degs_prop[,prop:=N/num_cells]
    degs_prop[,N_prop:=num_cells/sum(num_cells)]

    deg_prop_plot <- ggplot(data = degs_prop,
                        aes(x = factor(celltype), y = prop)) +
    geom_bar(stat = "identity") +
    labs(y = "Proportion of DEGs identified", x = "Cell Type", fill = "Cell Type") +
    geom_text(aes(label = scales::percent(prop, accuracy = 0.1)), vjust = -0.25, size = 3) +
    theme_cowplot() +
    theme(axis.text = element_text(size = 9)) +
    scale_fill_viridis(discrete = T)
    #save the graph to output_path
    suppressMessages(ggsave(path = output_path,
                            filename = "DEGs_proportion_per_cell_type.pdf",
                            plot=deg_prop_plot))
    suppressMessages(ggsave(path = output_path,
                            filename = "DEGs_proportion_per_cell_type.png",
                            plot=deg_prop_plot))


    #plot DEG direction and effect size
    celltype_DEGs_dt[,deg_direction:="Down"]
    celltype_DEGs_dt[logFC>0,deg_direction:="Up"]
    #count
    count_deg_direction_plot<-
        ggplot(data=celltype_DEGs_dt[,.N,by=.(celltype,deg_direction)],
               aes(x = factor(deg_direction), y = N, fill=factor(celltype)))+
        geom_bar(stat="identity") +
        facet_wrap(~factor(celltype))+
        ggtitle("Cell type counts of DEGs")+
        labs(y= "Number of DEGs", x = "DEG direction", fill="Cell Type") +
        theme_cowplot()+
        theme(axis.text = element_text(size=9))+
        scale_fill_viridis(discrete = T)
    #save the graph to output_path
    suppressMessages(
        ggsave(path = output_path,
               filename="DEGs_direction_count_per_cell_type.pdf",
               plot=count_deg_direction_plot))
    suppressMessages(
        ggsave(path = output_path,
               filename="DEGs_direction_count_per_cell_type.png",
               plot=count_deg_direction_plot))

    #plot volcano plots
    #get gene names - for sig genes
    genes <- unique(celltype_all_genes_dt[adj_pval<0.05,]$name)
    gene_IDs <- select(EnsDb.Hsapiens.v79, keys= genes, keytype = "GENEID", columns = c("GENEID","SYMBOL"))
    colnames(gene_IDs) <- c("ensembl_gene_id","hgnc_symbol")
    gene_IDs <- as.data.table(gene_IDs)
    setnames(gene_IDs,"ensembl_gene_id","name")
    #in case any ID's have more than one hgnc name found
    gene_IDs <- gene_IDs[!duplicated(gene_IDs$name),]
    setkey(gene_IDs,name)
    setkey(celltype_all_genes_dt,name)
    #append gene names
    celltype_all_genes_dt[, gene_name :=
                              gene_IDs[celltype_all_genes_dt, x.hgnc_symbol]]
    #add colour identifier
    celltype_all_genes_dt[,colour_ident:="Not Significant"]
    celltype_all_genes_dt[adj_pval<0.05 & logFC<0 ,
                            colour_ident:="Decreased DEG"]
    celltype_all_genes_dt[adj_pval<0.05 & logFC>0 ,
                            colour_ident:="Increased DEG"]

    cols <- c("Not Significant" = "grey", "Increased DEG" = "#FDE725FF",
              "Decreased DEG" = "#440154FF")
    #remove cell types without sig genes
    celltype_all_genes_dt <-
        celltype_all_genes_dt[(celltype %in%
                                    unique(celltype_DEGs_dt$celltype)),]

    volcano_plot<-
        ggplot(data=celltype_all_genes_dt,
               aes(x = logFC, y = -log10(adj_pval)))+ #scale FDR by log10
        geom_point(aes(colour=colour_ident),stat="identity",size=.2) +
        facet_wrap(~factor(celltype),scales = "free")+
        geom_hline(yintercept = -log10(0.05), colour="#990000",
                   linetype="dashed",size=.3) +
        #geom_vline(xintercept = 0, colour="black", size=.3) +
        labs(y= "-log10(FDR)", x = "Log2 Fold Change", colour="") +
        theme_cowplot()+
        theme(axis.text = element_text(size=9),axis.text.x=element_blank(),
              legend.text=element_text(size=10))+
        scale_colour_manual(values = cols)+
        geom_text_repel( #give top three genes per cell type by adjusted p value
            data = celltype_all_genes_dt[
                celltype_all_genes_dt[,.I[adj_pval %in%
                                              sort(adj_pval)[1:3]][1:3],
                                      by=celltype]$V1][!is.na(gene_name) & gene_name != ""],
            aes(label = gene_name),
            size = 3,
            box.padding = unit(0.35, "lines"),
            point.padding = unit(0.3, "lines")
        )
    #save the graph to output_path - increase size A4
    suppressMessages(ggsave(path = output_path,
                            filename = "volcano_plot_degs_cell_types.pdf",
                            plot=volcano_plot,
                            width = 210, height = 297, units = "mm"))
    suppressMessages(ggsave(path = output_path,
                            filename = "volcano_plot_degs_cell_types.png",
                            plot=volcano_plot,
                            width = 210, height = 297, units = "mm"))

    #plot boxplots of DEGs
    #check median log2 fold change for both up and down regulated genes
    #This will give a sense of directional effect size
    degs_median_lfc<-celltype_DEGs_dt[,median(logFC),by=deg_direction]
    #if only up or down reg DEGs add 0 for other
    if(nrow(degs_median_lfc)<2){
        dir <- c("Up","Down")
        degs_median_lfc<-
            rbind(degs_median_lfc,
                    data.table("deg_direction"=
                                   dir[dir!=degs_median_lfc$deg_direction],
                               "V1"=0))
        #set order so Up first
        setorder(degs_median_lfc,-deg_direction)
    }

    # Compute boxplot statistics
    calc_boxplot_stat <- function(x) {
        coef <- 1.5
        n <- sum(!is.na(x))
        # calculate quantiles
        stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
        names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
        iqr <- diff(stats[c(2, 4)])
        # set whiskers
        outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
        if (any(outliers)) {
            stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
        }
        return(stats)
    }

    degs_boxplot_plots<-
        ggplot(data=celltype_DEGs_dt,
               aes(x = factor(deg_direction), y = abs(logFC),
                    fill=factor(celltype)))+
        stat_summary(fun.data = calc_boxplot_stat, geom="boxplot") +
        facet_wrap(~factor(celltype),scales="free")+
        ggtitle("Cell type Log2 Fold Change of DEGs",
                subtitle = paste0("Median LFC of Down and Up regulated DEGs: ",
                                  round(degs_median_lfc$V1[[2]]*-1,2),", ",
                                  round(degs_median_lfc$V1[[1]],2)))+
        labs(y= "Log2 Fold Change", x = "DEG direction", fill="Cell Type") +
        theme_cowplot()+
        theme( axis.text = element_text(size=9))+
        scale_fill_viridis(discrete = T)
    #save the graph to output_path
    suppressMessages(ggsave(path = output_path,
                            filename = "deg_boxplots_cell_types.pdf",
                            plot=degs_boxplot_plots))
    suppressMessages(ggsave(path = output_path,
                            filename = "deg_boxplots_cell_types.png",
                            plot=degs_boxplot_plots))
}
