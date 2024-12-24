#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)
library(Seurat)
library(MAST)

#run diff expr
seurat_obj <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Seurat_Objects/8_processed_seurat_obj_integrated_azimuthAnnotations_ourAnnotations_doubletsFiltered_FCD1Filtered.RDS")
tle_control_samples <- c("TLE36", "TLE6", "TLE7", "FC10101")
seurat_obj$status_tle <- seurat_obj$status
seurat_obj$status_tle[which(seurat_obj$cleaned_ident %in% tle_control_samples)] <- "tle_control"

DEG_res <- FindMarkers(seurat_obj, ident.1 = "case", ident.2 = "control", group.by = "status_tle", min.pct = 0.10, logfc.threshold = 0, max.cells.per.ident = 30000)
saveRDS(DEG_res, paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Pseudobulk/wilcox_10/All/DEG_res.RDS"))
# DEG_res <- FindMarkers(seurat_obj, ident.1 = "case", ident.2 = "control", group.by = "status_tle", test.use = "MAST", latent.vars = c("orig.ident"), min.pct = 0.10, logfc.threshold = 0, max.cells.per.ident = 30000)
# saveRDS(DEG_res, paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Pseudobulk/mast_10/All/DEG_res.RDS"))
