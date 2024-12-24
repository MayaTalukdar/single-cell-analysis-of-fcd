#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)
library(Seurat)
library(MAST)

cell_type = commandArgs(trailingOnly=TRUE)[1]
cleaned_cell_types <- c("L5/6_NP" = "L5/6 NP", "L5_ET" = "L5 ET", "L6_IT_Car3" = "L6 IT Car3", "L6_CT" = "L6 CT", "L6_IT" = "L6 IT", "L2/3_IT" = "L2/3 IT", "L4_IT" = "L4 IT", "L5_IT" = "L5 IT")
if (cell_type %in% names(cleaned_cell_types))
{
    cell_type <- cleaned_cell_types[cell_type]
}
print(cell_type)

#run diff expr
seurat_obj <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Seurat_Objects/9_processed_seurat_obj_integrated_azimuthAnnotations_ourAnnotations_doubletsFiltered_FCD1Filtered_looseGenotyping.RDS")
Idents(seurat_obj) <- "final_annots"

DEG_res <- FindMarkers(seurat_obj, ident.1 = "hetero", ident.2 = "ref-hom", group.by = "genotype", subset.ident = cell_type, min.pct = 0.10, logfc.threshold = 0)
saveRDS(DEG_res, paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Het_Ref/wilcox_10/All/", gsub("/", "_", gsub(" ", "_", cell_type)), "_DEG_res.RDS"))
DEG_res <- FindMarkers(seurat_obj, ident.1 = "hetero", ident.2 = "ref-hom", group.by = "genotype", subset.ident = cell_type, test.use = "MAST", latent.vars = c("orig.ident"), min.pct = 0.10, logfc.threshold = 0)
saveRDS(DEG_res, paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Het_Ref/mast_10/All/", gsub("/", "_", gsub(" ", "_", cell_type)), "_DEG_res.RDS"))
