library(tidyverse)
library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")
# library(msigdbr)
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(DOSE)
# library(enrichplot)
# library(limma)
# library(gplots)
# library(marray)
# library(RMySQL)
# library(stringr)
# library(reshape2)
# library(dplyr)
# library(fgsea)
# library(pheatmap)
# library(cowplot)
# library(patchwork)
# library(scCustomize)
# library(fgsea)
# library(lsa)
# library(rjson)
# library(harmony)
# library(Azimuth)

#source("/n/groups/walsh/indData/Maya/microCHIP_AD_project/GOT/FCD_Analysis/Final_Analysis_For_First_Submission/util.R")

# #######################
# #CREATE UNPROCESSED SEURAT OBJECT
# #######################
# data.dirs <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/FCD_project/raw_data/paths_to_fcd_files_cellbender_output.txt", header = FALSE)))

# #read in all 10X files as matrices 
# list.of.10X.files <- list()
# for (dir in data.dirs)
# {
#   mat <-  Read10X(data.dir=paste0(dir, "/outs/cellbender_feature_bc_matrix/"))
#   sample_name <- gsub("_intron_count", "", basename(dir))
#   print(sample_name)
#   colnames(mat) <- paste0(sample_name, "-", colnames(mat))
#   list.of.10X.files[[sample_name]] <- mat
# }
# print("READ IN ALL MATRICES!")

# #create seurat objects from matrices 
# list.of.Seurat.objs <- lapply(seq_along(list.of.10X.files), function(x) CreateSeuratObject(counts = list.of.10X.files[[x]], min.cells = 0, min.features = 10, project = names(list.of.10X.files)[x])) 
# names(list.of.Seurat.objs) <- names(list.of.10X.files)
# print("CREATED ALL SEURAT OBJECTS!")

# #create merged seurat object 
# merged.Seurat.obj <- merge(list.of.Seurat.objs[[1]], y = list.of.Seurat.objs[2:length(list.of.Seurat.objs)])
# print("MERGED SEURAT OBJECT!")

# #save merged seurat object 
# saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/Seurat_Objects/0_unprocessed_Seurat_obj.RDS")
# print("SAVE OBJECT!")

# ######################
# #DETERMINE DOUBLETS WITH SCRUBLET
# ######################
# merged.Seurat.obj <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/Seurat_Objects/0_unprocessed_Seurat_obj.RDS")
# prop_doublets <- list()
# num_doublets <- list()
# num_singlets <- list()
# cells_to_keep <- c()

# #chung 
# unique_idents <- list.files("/n/groups/walsh/indData/Maya/FCD_project/raw_data/Chung/AlignedFiles/", pattern = "SRR")
# unique_idents <- setdiff(unique_idents, "SRR24190238")
# for (ident in unique_idents)
# {
#   #set up sample names 
#   sample_name <- gsub("_intron_count", "", ident)
#   print(sample_name)

#   #read in scrublet scores
#   scrublet_scores <- unname(unlist(read.csv(paste0("/n/groups/walsh/indData/Maya/FCD_project/raw_data/Chung/Scrublet/Output/", ident, "_scrublet_predictions.csv"), header = TRUE)))
#   names(scrublet_scores) <- unname(unlist(read.table(paste0("/n/groups/walsh/indData/Maya/FCD_project/raw_data/Chung/AlignedFiles/", ident, "/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"))))
#   prop_doublets[[sample_name]] <- length(which(scrublet_scores == "True"))/length(scrublet_scores)
#   num_doublets[[sample_name]] <- length(which(scrublet_scores == "True"))
#   num_singlets[[sample_name]] <- length(which(scrublet_scores == "False"))

#   #get cells to keep 
#   current_cells_to_keep <- sapply(names(scrublet_scores)[which(scrublet_scores == "False")], function(x) paste0(sample_name, "-", x))
#   current_cells_to_keep <- unname(gsub("_", "-", current_cells_to_keep))
#   cells_to_keep <- c(cells_to_keep, current_cells_to_keep)
# }

# #ours
# wds <- unlist(unname(read.table("/n/groups/walsh/indData/Maya/FCD_project/raw_data/paths_to_fcd_files_cellbender_output.txt", header = FALSE)))
# wds <- wds[-which(grepl("SRR", wds))]
# for (wd in wds)
# {
#     #set up sample name 
#     ident <- basename(wd)
#     sample_name <- gsub("_intron_count", "", ident)
#     print(sample_name)

#     #read in scrublet scores 
#     scrublet_scores <- unname(unlist(read.csv(paste0("/n/groups/walsh/indData/Maya/FCD_project/raw_data/Scrublet/Output/", ident, "_scrublet_predictions.csv"), header = TRUE)))
#     names(scrublet_scores) <- unname(unlist(read.table(paste0(wd, "/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"))))
#     prop_doublets[[sample_name]] <- length(which(scrublet_scores == "True"))/length(scrublet_scores)
#     num_doublets[[sample_name]] <- length(which(scrublet_scores == "True"))
#     num_singlets[[sample_name]] <- length(which(scrublet_scores == "False"))

#     #get cells to keep 
#     current_cells_to_keep <- sapply(names(scrublet_scores)[which(scrublet_scores == "False")], function(x) paste0(sample_name, "-", x))
#     current_cells_to_keep <- unname(gsub("_", "-", current_cells_to_keep))
#     cells_to_keep <- c(cells_to_keep, current_cells_to_keep)
#  }

# #make sure that we have cells to retain from all samples 
# length(intersect(cells_to_keep, colnames(merged.Seurat.obj)))/ncol(merged.Seurat.obj) #~83% cells are retained

# #get all putative doublets
# cells_to_keep <- intersect(unname(cells_to_keep), colnames(merged.Seurat.obj))
# saveRDS(cells_to_keep, "/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/Seurat_Objects/cells_not_identified_as_doublets_by_scrublet.RDS")

# #add this metadata
# merged.Seurat.obj$is_putative_doublet <- (!row.names(merged.Seurat.obj@meta.data) %in% cells_to_keep)
# saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/Seurat_Objects/1_processed_seurat_obj_with_scrublet_info.RDS")

# ###################
# #QUALITY CONTROL
# ###################
# merged.Seurat.obj <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/Seurat_Objects/1_processed_seurat_obj_with_scrublet_info.RDS")
# unique(merged.Seurat.obj$orig.ident)
# alt_ids <- read.csv("/n/groups/walsh/indData/Maya/FCD_project/raw_data/alt_ids.csv", header = FALSE)
# alt_ids <- setNames(alt_ids$V2, alt_ids$V1)
# merged.Seurat.obj$cleaned_ident <- merged.Seurat.obj$orig.ident
# merged.Seurat.obj$cleaned_ident[which(merged.Seurat.obj$cleaned_ident %in% names(alt_ids))] <- unname(alt_ids[merged.Seurat.obj$cleaned_ident[which(merged.Seurat.obj$cleaned_ident %in% names(alt_ids))]])

# #add metrics using scCustomize
# merged.Seurat.obj<- Add_Mito_Ribo_Seurat(seurat_object = merged.Seurat.obj, species = "Human")
# merged.Seurat.obj<- Add_Cell_Complexity_Seurat(seurat_object = merged.Seurat.obj)
# all.equal(colnames(merged.Seurat.obj), row.names(merged.Seurat.obj@meta.data))

# #create plots of qc metrics before filtering
# pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/Quality_Control/qc_plots_pre_filtering_with_thresholds.pdf", width = 13, height = 15)
# Idents(merged.Seurat.obj) <- "cleaned_ident"
# p1 <- QC_Plots_Genes(seurat_object = merged.Seurat.obj, plot_title = "Genes Per Cell", low_cutoff = 300, high_cutoff = 8000, raster = TRUE)
# p2 <- QC_Plots_UMIs(seurat_object = merged.Seurat.obj, plot_title = "UMIs Per Cell", low_cutoff = 300, high_cutoff = 25000, raster = TRUE)
# p3 <- QC_Plots_Mito(seurat_object = merged.Seurat.obj, plot_title = "Mito Gene % Per Cell", high_cutoff = 5, raster = TRUE)
# p4 <- QC_Plots_Feature(seurat_object = merged.Seurat.obj, feature = "percent_ribo", plot_title = "Ribo Gene % Per Cell", high_cutoff = 5, raster = TRUE)
# wrap_plots(p1, p2, p3, p4, ncol = 1)
# dev.off()

# #determine cells that will be lost by each qc filtering step 
# cells_filtered_pct_mito <- row.names(merged.Seurat.obj@meta.data %>% filter(percent_mito >= 5))
# cells_filtered_pct_ribo <- row.names(merged.Seurat.obj@meta.data %>% filter(percent_ribo >= 5))
# cells_filtered_nFeatures <- row.names(merged.Seurat.obj@meta.data %>% filter(nFeature_RNA <= 300))
# cells_filtered_nCounts <- row.names(merged.Seurat.obj@meta.data %>% filter(nCount_RNA <=300))
# filtering_list <- list("pct_mito" = cells_filtered_pct_mito, "pct_ribo" = cells_filtered_pct_ribo, "nFeatures" = cells_filtered_nFeatures, "nCounts" = cells_filtered_nCounts)
# print("NUMBER OF CELLS FILTERED PER CRITERIA: ")
# print(lapply(filtering_list, length))
# print(length(unique(do.call(c, filtering_list))))
# print(length(unique(do.call(c, filtering_list)))/ncol(merged.Seurat.obj))

# #filter cells 
# merged.Seurat.obj<- subset(merged.Seurat.obj, subset = percent_mito < 5)
# merged.Seurat.obj <- subset(merged.Seurat.obj, subset = percent_ribo < 5)
# merged.Seurat.obj <- subset(merged.Seurat.obj, subset = nFeature_RNA > 300)
# merged.Seurat.obj <- subset(merged.Seurat.obj, subset = nCount_RNA > 300)
# saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/Seurat_Objects/2_processed_seurat_obj_filtered.RDS")

# #####################################
# #NORMALIZATION & DIMENSIONALITY REDUCTION
# ##normalize and residualize cells 
# #####################################
# merged.Seurat.obj <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/Seurat_Objects/2_processed_seurat_obj_filtered.RDS")

# #convert to seurat assay v5
# merged.Seurat.obj[["RNA5"]] <- as(object = merged.Seurat.obj[["RNA"]], Class = "Assay5")

# #normalize data
# DefaultAssay(merged.Seurat.obj) <- "RNA5"
# Idents(merged.Seurat.obj) <- "cleaned_ident"
# merged.Seurat.obj[["RNA5"]] <- split(merged.Seurat.obj[["RNA5"]], f = merged.Seurat.obj$orig.ident)
# merged.Seurat.obj<- NormalizeData(merged.Seurat.obj)
# merged.Seurat.obj <- FindVariableFeatures(merged.Seurat.obj, selection.method = "vst", nfeatures = 3000)
# merged.Seurat.obj <- ScaleData(merged.Seurat.obj, vars.to.regress = c("percent_mito", "percent_ribo", "nCount_RNA", "nFeature_RNA"))

# #perform dimensionality reduction
# merged.Seurat.obj <- RunPCA(merged.Seurat.obj)
# merged.Seurat.obj <- ProjectDim(object = merged.Seurat.obj)
# merged.Seurat.obj <- FindNeighbors(merged.Seurat.obj, dims = 1:30, reduction = "pca")
# merged.Seurat.obj <- FindClusters(merged.Seurat.obj, resolution = 1.2, cluster.name = "unintegrated_clusters.1.2")
# merged.Seurat.obj <- FindClusters(merged.Seurat.obj, resolution = 0.6, cluster.name = "unintegrated_clusters.0.6")
# merged.Seurat.obj <- RunUMAP(merged.Seurat.obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# #look at unintegrated cell types
# DefaultAssay(merged.Seurat.obj) <- "RNA5"
# pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/UMAPs/unintegrated_umap.pdf", width = 12)
# DimPlot(merged.Seurat.obj, reduction = "umap.unintegrated", group.by = c("orig.ident", "cleaned_ident"), label = TRUE) + NoLegend()
# dev.off()

# #save
# saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/Seurat_Objects/3_processed_seurat_obj_normalized_and_dimRed.RDS")

# #####################################
# #INTEGRATE
# #####################################
# #read in objects
# merged.Seurat.obj <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/Seurat_Objects/3_processed_seurat_obj_normalized_and_dimRed.RDS")
# DefaultAssay(merged.Seurat.obj) <- "RNA5"

# merged.Seurat.obj  <- IntegrateLayers(
#   object = merged.Seurat.obj , method = RPCAIntegration,
#   orig.reduction = "pca", new.reduction = "integrated.rpca",
#   verbose = TRUE
# )
# print("Completedd rpca!")
# saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/Seurat_Objects/4_processed_seurat_obj_integrated_rpca.RDS")

# # integrate 
# merged.Seurat.obj  <- IntegrateLayers(
#   object = merged.Seurat.obj , method = CCAIntegration,
#   orig.reduction = "pca", new.reduction = "integrated.cca",
#   verbose = TRUE
# )
# print("Completed cca!")
# saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/Seurat_Objects/4_processed_seurat_obj_integrated_rpca_cca.RDS")

# merged.Seurat.obj <- IntegrateLayers(
#  object = merged.Seurat.obj, method = HarmonyIntegration,
#   orig.reduction = "pca", new.reduction = "harmony",
#   verbose = TRUE
#  )
# print("Completed harmony!")
# saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/Seurat_Objects/4_processed_seurat_obj_integrated_rpca_cca_harmony.RDS")

# #####################################
# #ANNOTATION
# #####################################
# # merged.Seurat.obj <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/Seurat_Objects/4_processed_seurat_obj_integrated_rpca_cca_harmony.RDS")

# # #run Azimuth
# # merged.Seurat.obj <- JoinLayers(merged.Seurat.obj)
# # merged.Seurat.obj <- RunAzimuth(merged.Seurat.obj, reference = "humancortexref")
# # saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/Seurat_Objects/5_processed_seurat_obj_integrated_azimuthAnnotations.RDS")

# #####################################
# #CLUSTERING ON INTEGRATED REDUCTIONS
# #####################################
# merged.Seurat.obj <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/Seurat_Objects/5_processed_seurat_obj_integrated_azimuthAnnotations.RDS")
# merged.Seurat.obj <- JoinLayers(merged.Seurat.obj)


# #cca
# merged.Seurat.obj <- FindNeighbors(merged.Seurat.obj, reduction = "integrated.cca", dims = 1:30)
# merged.Seurat.obj <- FindClusters(merged.Seurat.obj, resolution = 1, cluster.name = "cca_clusters")
# merged.Seurat.obj <- RunUMAP(merged.Seurat.obj, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
# # umap_by_cleaned_ident <- DimPlot(merged.Seurat.obj, "umap.cca", group.by = "cleaned_ident", label = FALSE)
# # umap_by_predicted_class <- DimPlot(merged.Seurat.obj, "umap.cca", group.by = "predicted.class", label = TRUE, label.size = 2) + NoLegend()
# # umap_by_predicted_subclass <- DimPlot(merged.Seurat.obj, "umap.cca", group.by = "predicted.subclass", label = TRUE, label.size = 2) + NoLegend()
# # umap_by_doublet <- DimPlot(merged.Seurat.obj, "umap.cca", group.by = "is_putative_doublet", label = FALSE) 

# # pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/UMAPs/integrated_umaps_cca.pdf", width = 24, height = 6)
# # wrap_plots(umap_by_cleaned_ident, umap_by_predicted_class, umap_by_predicted_subclass, umap_by_doublet)
# # dev.off()

# #rpca
# merged.Seurat.obj <- FindNeighbors(merged.Seurat.obj, reduction = "integrated.rpca", dims = 1:30)
# merged.Seurat.obj <- FindClusters(merged.Seurat.obj, resolution = 1, cluster.name = "rpca_clusters")
# merged.Seurat.obj <- RunUMAP(merged.Seurat.obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
# # umap_by_cleaned_ident <- DimPlot(merged.Seurat.obj, "umap.rpca", group.by = "cleaned_ident", label = FALSE)
# # umap_by_predicted_class <- DimPlot(merged.Seurat.obj, "umap.rpca", group.by = "predicted.class", label = TRUE, label.size = 2) + NoLegend()
# # umap_by_predicted_subclass <- DimPlot(merged.Seurat.obj, "umap.rpca", group.by = "predicted.subclass", label = TRUE, label.size = 2) + NoLegend()
# # umap_by_doublet <- DimPlot(merged.Seurat.obj, "umap.rpca", group.by = "is_putative_doublet", label = FALSE) 

# # pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/UMAPs/integrated_umaps_rpca.pdf", width = 24, height = 6)
# # wrap_plots(umap_by_cleaned_ident, umap_by_predicted_class, umap_by_predicted_subclass, umap_by_doublet)
# # dev.off()

# #harmony
# merged.Seurat.obj <- FindNeighbors(merged.Seurat.obj, reduction = "harmony", dims = 1:30)
# merged.Seurat.obj <- FindClusters(merged.Seurat.obj, resolution = 1, cluster.name = "harmony_clusters")
# merged.Seurat.obj <- RunUMAP(merged.Seurat.obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
# # umap_by_cleaned_ident <- DimPlot(merged.Seurat.obj, "umap.harmony", group.by = "cleaned_ident", label = FALSE)
# # umap_by_predicted_class <- DimPlot(merged.Seurat.obj, "umap.harmony", group.by = "predicted.class", label = TRUE, label.size = 2) + NoLegend()
# # umap_by_predicted_subclass <- DimPlot(merged.Seurat.obj, "umap.harmony", group.by = "predicted.subclass", label = TRUE, label.size = 2) + NoLegend()
# # umap_by_doublet <- DimPlot(merged.Seurat.obj, "umap.harmony", group.by = "is_putative_doublet", label = FALSE) 

# # pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/UMAPs/integrated_umaps_harmony.pdf", width = 24, height = 6)
# # wrap_plots(umap_by_cleaned_ident, umap_by_predicted_class, umap_by_predicted_subclass, umap_by_doublet)
# # dev.off()

# saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/Seurat_Objects/6_processed_seurat_obj_integrated_withClustering.RDS")

#####################################
#CLUSTERING ON INTEGRATED REDUCTIONS
#####################################
seurat_obj <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/Seurat_Objects/6_processed_seurat_obj_integrated_withClustering.RDS")

#filter out doublts
Idents(seurat_obj) <- "is_putative_doublet"
seurat_obj <- subset(seurat_obj, idents = FALSE)
table(seurat_obj$is_putative_doublet)

#filter out Chung FCD1 sample
Idents(seurat_obj) <- "cleaned_ident"
seurat_obj <- subset(seurat_obj, idents = "SRR24190237", invert = TRUE)

saveRDS(seurat_obj, "/n/groups/walsh/indData/Maya/FCD_project/analysis/1_Analyze_FCD_Object/Output/Seurat_Objects/7_processed_seurat_obj_integrated_withClustering_doubletsFiltered_FCD1Filtered.RDS")
