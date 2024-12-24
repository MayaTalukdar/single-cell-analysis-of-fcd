#!/usr/bin/env Rscript
###############################
#I/O
###############################
library(tidyverse)
library(ggplot2)
library(Seurat)
library(MAST)

#add in info about tle status
seurat_obj <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Seurat_Objects/8_processed_seurat_obj_integrated_azimuthAnnotations_ourAnnotations_doubletsFiltered_FCD1Filtered.RDS")
tle_control_samples <- c("TLE36", "TLE6", "TLE7", "FC10101")
seurat_obj$status_tle <- seurat_obj$status
seurat_obj$status_tle[which(seurat_obj$cleaned_ident %in% tle_control_samples)] <- "tle_control"

#add in info about cell groupings
#add in higher level annots
higher_level_annots <- c(
"Oligo" = "Glial", 
"OPC" = "Glial", 
"Astro" = "Glial", 
"L2/3 IT" = "Neuron",
"Micro-PVM" = "Glial", 
"Pvalb" = "Neuron", 
"L4 IT" = "Neuron", 
"L5 IT" = "Neuron", 
"L5 ET" = "Neuron", 
"Vip" = "Neuron", 
"Sst" = "Neuron", 
"Sncg" = "Neuron",
"Lamp5" = "Neuron", 
"L5/6 NP" = "Neuron", 
"L6 CT" = "Neuron", 
"L6 IT" = "Neuron", 
"L6 IT Car3" = "Neuron", 
"L6b" = "Neuron", 
"Endo" = "Glial", 
"VLMC" = "Glial"
)
seurat_obj$final_class <- setNames(higher_level_annots[seurat_obj$final_annots], names(seurat_obj$final_annots))
Idents(seurat_obj) <- "final_class"

DEG_res <- FindMarkers(seurat_obj, ident.1 = "case", ident.2 = "control", group.by = "status_tle", subset.ident = "Neuron",  min.pct = 0.10, logfc.threshold = 0, max.cells.per.ident = 30000)
saveRDS(DEG_res, paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_NonTLEControl_Neurons/wilcox_10/All/DEG_res.RDS"))
# DEG_res <- FindMarkers(seurat_obj, ident.1 = "case", ident.2 = "control", group.by = "status_tle", subset.ident = "Neuron", test.use = "MAST", latent.vars = c("orig.ident"), min.pct = 0.10, logfc.threshold = 0, max.cells.per.ident = 30000)
# saveRDS(DEG_res, paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_NonTLEControl_Neurons/mast_10/All/DEG_res.RDS"))
