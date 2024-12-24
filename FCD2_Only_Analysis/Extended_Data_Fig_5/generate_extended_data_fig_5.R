library(tidyverse)
library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(limma)
library(gplots)
library(marray)
library(RMySQL)
library(stringr)
library(reshape2)
library(dplyr)
library(fgsea)
library(pheatmap)
library(cowplot)
library(patchwork)
library(scCustomize)
library(fgsea)
library(lsa)
library(rjson)
library(harmony)
library(Azimuth)

source("/n/groups/walsh/indData/Maya/microCHIP_AD_project/GOT/FCD_Analysis/Final_Analysis_For_First_Submission/util.R")

#################
#PROCESS INDIVIDUAL SAMPLES
#################
samples <- c( "e348", "e174", "e254", "e274", "e174")
for (sample in samples)
{
    print("****************")
    print(sample)
    print("****************")

    output_dir <- paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/Per_Sample_Analysis/", sample)
    dir.create(output_dir)

    #######################
    #CREATE UNPROCESSED SEURAT OBJECT
    #######################
    #read in data 
    data.dirs <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/FCD_project/raw_data/paths_to_fcd_files_cellbender_output.txt", header = FALSE)))
    dir <- data.dirs[which(grepl(sample, data.dirs))]
    mat <-  Read10X(data.dir=paste0(dir, "/outs/cellbender_feature_bc_matrix/"))
    sample_name <- gsub("_intron_count", "", basename(dir))
    print(sample_name)
    colnames(mat) <- paste0(sample_name, "-", colnames(mat))
    seurat_obj <- CreateSeuratObject(counts =  mat, min.cells = 0, min.features = 10, project = sample)
    seurat_obj$nCount_RNA <- colSums(seurat_obj)
    seurat_obj$nFeature_RNA <- colSums(x = GetAssayData(object = seurat_obj, slot = "counts") > 0)  # nFeatureRNA

    #add in metadata from the full object around doublets and then filter them 
    putative_singlets <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/Seurat_Objects/cells_not_identified_as_doublets_by_scrublet.RDS")
    seurat_obj <- seurat_obj[,which(colnames(seurat_obj) %in% putative_singlets)]

    #add in cleaned ident
    alt_ids <- read.csv("/n/groups/walsh/indData/Maya/FCD_project/raw_data/alt_ids.csv", header = FALSE)
    alt_ids <- setNames(alt_ids$V2, alt_ids$V1)
    seurat_obj$cleaned_ident <- seurat_obj$orig.ident
    seurat_obj$cleaned_ident[which(seurat_obj$cleaned_ident %in% names(alt_ids))] <- unname(alt_ids[seurat_obj$cleaned_ident[which(seurat_obj$cleaned_ident %in% names(alt_ids))]])
    
    ###################
    #QUALITY CONTROL
    ###################
    qc_dir <- paste0(output_dir, "/Quality_Control/")
    dir.create(qc_dir)
    
    #add metrics using scCustomize
    seurat_obj<- Add_Mito_Ribo_Seurat(seurat_object = seurat_obj, species = "Human")
    all.equal(colnames(seurat_obj), row.names(seurat_obj@meta.data))

    #create plots of qc metrics before filtering
    pdf(paste0(qc_dir, "qc_plots_pre_filtering_with_thresholds.pdf"), width = 13, height = 15)
    Idents(seurat_obj) <- "cleaned_ident"
    p1 <- QC_Plots_Genes(seurat_object = seurat_obj, plot_title = "Genes Per Cell", low_cutoff = 300, high_cutoff = 8000, raster = TRUE)
    p2 <- QC_Plots_UMIs(seurat_object = seurat_obj, plot_title = "UMIs Per Cell", low_cutoff = 300, high_cutoff = 25000, raster = TRUE)
    p3 <- QC_Plots_Mito(seurat_object = seurat_obj, plot_title = "Mito Gene % Per Cell", high_cutoff = 5, raster = TRUE)
    p4 <- QC_Plots_Feature(seurat_object = seurat_obj, feature = "percent_ribo", plot_title = "Ribo Gene % Per Cell", high_cutoff = 5, raster = TRUE)
    print(wrap_plots(p1, p2, p3, p4, ncol = 1))
    dev.off()

    #determine cells that will be lost by each qc filtering step 
    cells_filtered_pct_mito <- row.names(seurat_obj@meta.data %>% filter(percent_mito >= 5))
    cells_filtered_pct_ribo <- row.names(seurat_obj@meta.data %>% filter(percent_ribo >= 5))
    cells_filtered_nFeatures <- row.names(seurat_obj@meta.data %>% filter(nFeature_RNA <= 300))
    cells_filtered_nCounts <- row.names(seurat_obj@meta.data %>% filter(nCount_RNA <=300))
    filtering_list <- list("pct_mito" = cells_filtered_pct_mito, "pct_ribo" = cells_filtered_pct_ribo, "nFeatures" = cells_filtered_nFeatures, "nCounts" = cells_filtered_nCounts)
    print("NUMBER OF CELLS FILTERED PER CRITERIA: ")
    print(lapply(filtering_list, length))
    print(length(unique(do.call(c, filtering_list))))
    print(length(unique(do.call(c, filtering_list)))/ncol(seurat_obj))

    #filter cells 
    seurat_obj<- subset(seurat_obj, subset = percent_mito < 5)
    seurat_obj <- subset(seurat_obj, subset = percent_ribo < 5)
    seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 300)
    seurat_obj <- subset(seurat_obj, subset = nCount_RNA > 300)

    #####################################
    #NORMALIZATION & DIMENSIONALITY REDUCTION
    ##normalize and residualize cells 
    #####################################
    umap_dir <-  paste0(output_dir, "/UMAPs/")
    dir.create(umap_dir)

    #convert to seurat assay v5
    seurat_obj[["RNA5"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay5")

    #normalize data
    DefaultAssay(seurat_obj) <- "RNA5"
    Idents(seurat_obj) <- "cleaned_ident"
    seurat_obj[["RNA5"]] <- split(seurat_obj[["RNA5"]], f = seurat_obj$orig.ident)
    seurat_obj<- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
    seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("percent_mito", "percent_ribo", "nCount_RNA", "nFeature_RNA"))

    #perform dimensionality reduction
    seurat_obj <- RunPCA(seurat_obj)
    seurat_obj <- ProjectDim(object = seurat_obj)
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, reduction = "pca")
    seurat_obj <- FindClusters(seurat_obj, resolution = 1.2, cluster.name = "unintegrated_clusters.1.2")
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.6, cluster.name = "unintegrated_clusters.0.6")
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

    #look at unintegrated cell types
    DefaultAssay(seurat_obj) <- "RNA5"
    pdf(paste0(umap_dir, "unintegrated_umap.pdf"), width = 12)
    print(DimPlot(seurat_obj, reduction = "umap.unintegrated", group.by = c("orig.ident", "cleaned_ident"), label = TRUE) + NoLegend())
    dev.off()

    #run Azimuth
    seurat_obj <- RunAzimuth(seurat_obj, reference = "humancortexref")
    pdf(paste0(umap_dir, "azimuth_annotations.pdf"), width = 12)
    p1 <- DimPlot(seurat_obj, reduction = "umap.unintegrated", group.by = c("unintegrated_clusters.0.6"), label = TRUE, label.size = 2, repel = FALSE) + NoLegend()
    p2 <- DimPlot(seurat_obj, reduction = "umap.unintegrated", group.by = c("predicted.subclass"), label = TRUE, label.size = 2, repel = TRUE) + NoLegend()
    print(wrap_plots(p1, p2, nrow = 1))
    dev.off()

    #create violin plots of prediction scores
    label_transfer_score_df <- seurat_obj@meta.data %>% 
    group_by(unintegrated_clusters.0.6, predicted.subclass) %>% 
    summarize(med_pred_score = median(predicted.subclass.score), .groups = "drop") %>% 
    pivot_wider(names_from = predicted.subclass, 
                values_from = med_pred_score, 
                values_fill = 0) %>% as.data.frame()
    row.names(label_transfer_score_df) <- paste0("Cluster ", label_transfer_score_df$unintegrated_clusters.0.6)
    label_transfer_score_df <- label_transfer_score_df[,-1]

    label_transfer_long <- label_transfer_score_df %>%
    rownames_to_column(var = "Cluster") %>%
    pivot_longer(cols = -Cluster, 
                names_to = "Cell_Type", 
                values_to = "Score")

    label_transfer_long$Cell_Type <- factor(label_transfer_long$Cell_Type, 
                                                levels = sort(unique(label_transfer_long$Cell_Type)))

    max_cell_type <- sapply(unique(label_transfer_long$Cluster), function(x) max((label_transfer_long %>% filter(Cluster == x))$Score))
    max_cell_type <- max_cell_type[rev(order(max_cell_type))]
    label_transfer_long$Cluster <- factor(label_transfer_long$Cluster, levels = names(max_cell_type))

    pdf(paste0(umap_dir, "predicted_subclass_vlnplot.pdf"), width = 10, height = 10)
    print(ggplot(label_transfer_long, aes(x = Cluster, y = Score)) +
    geom_violin(fill = "gray80", alpha = 0.7, trim = TRUE) + # Uniform gray violin fill
    geom_jitter(aes(color = Cell_Type), width = 0.2, size = 1, alpha = 0.7) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Cluster", y = "Prediction Score") +
    geom_hline(yintercept = 0.6, color = "red", linetype = "dashed"))
    dev.off()
}





