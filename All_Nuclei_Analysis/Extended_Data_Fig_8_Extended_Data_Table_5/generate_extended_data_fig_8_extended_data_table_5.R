###########################
#I/O
###########################
library(tidyverse)
library(ggplot2)
library(Seurat)
library(fgsea)
library(pheatmap)
library(cowplot)
library(patchwork)
library(scCustomize)
library(CellChat)
library(msigdbr)
library(org.Hs.eg.db)
library(readxl)
library(openxlsx)
library(ggrepel)
library(gridExtra)

mat_for_plot_list <- list()

#####################
#CASE/CONTROL GSEA (GLIAL)
#####################
#write out version of DEG files with only protein coding genes
protein_coding_genes <- setNames(unlist(data.table::fread("/n/groups/walsh/indData/Maya/sex_differences_in_the_heart/human_protein_coding_genes.txt", header = TRUE)), NULL)
dir.create("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_NonTLEControl_Glial/wilcox_10/Protein-Coding/")
for (file in list.files("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_NonTLEControl_Glial/wilcox_10/All/", pattern = ".RDS"))
{
	res <- readRDS(paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_NonTLEControl_Glial/wilcox_10/All/", file))
	res <- res %>% filter(row.names(res) %in% protein_coding_genes)
	saveRDS(res, paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_NonTLEControl_Glial/wilcox_10/Protein-Coding/", file))
  write.csv(res, paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_NonTLEControl_Glial/wilcox_10/Protein-Coding/", gsub(".RDS", "", file), ".csv"))
}

#get hallmark gene sets
h_gene_sets = msigdbr(species = "human", category = "H")
hallmark_gs <- lapply(unique(h_gene_sets$gs_name), function(x) h_gene_sets$gene_symbol[which(h_gene_sets$gs_name == x)])
names(hallmark_gs) <- unique(h_gene_sets$gs_name)
pathway_list <- hallmark_gs

#add epilepsy gene sets
##DisGeNET
epilepsy_gene_set <- jsonlite::fromJSON(txt="/n/groups/walsh/indData/Maya/FCD_project/raw_data/Ontologies/disgenet_epilepsy_json_string.json")
epilepsy_gene_set <-  epilepsy_gene_set$associations$gene$symbol
pathway_list[["epilepsy_DisGeNET"]] <- epilepsy_gene_set

##Macnee et al., 
macnee_data <- read.csv("/n/groups/walsh/indData/Maya/FCD_project/raw_data/Ontologies/macnee_epilepsy_genes.csv", header = TRUE)
pathway_list[["epilepsy_macnee_high_conf"]] <- (macnee_data %>% filter(Classification == "tier 1"))$Symbol
pathway_list[["epilepsy_macnee"]] <- macnee_data$Symbol

#add synapse gene sets 
##SynGO
syngo_list <- read_excel("/n/groups/walsh/indData/Maya/FCD_project/raw_data/Ontologies/SynGo/syngo_annotations.xlsx") %>% as.data.frame()
syngo_list <- split(syngo_list, syngo_list$go_name)
syngo_list <- lapply(syngo_list, function(x) x$hgnc_symbol)
names(syngo_list) <- sapply(names(syngo_list), function(x) strsplit(x, "(", fixed = TRUE)[[1]][1])
pathway_list <- c(pathway_list , syngo_list) 

#add autism gene set 
sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv", header = TRUE)
SFARI_genes <- unique(sfari$gene.symbol)
high_conf_SFARI_genes <- unique((sfari %>% filter(gene.score == 1))$gene.symbol)
pathway_list[["sfari"]] <- SFARI_genes
pathway_list[["high_conf_sfari"]] <- high_conf_SFARI_genes

#run gsea
wd <- "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_NonTLEControl_Glial/wilcox_10/Protein-Coding"
gsea_res_list <- list()
for (file in list.files(wd, ".RDS"))
{
  type <- gsub("_DEG_res.RDS", "", file)
  print(type)
  DEG_res <- readRDS(paste0(wd, "/", file)) %>% arrange(desc(avg_log2FC)) %>% filter(pct.1 > 0.10) %>% filter(pct.2 > 0.10)
  ranks <- setNames(DEG_res$avg_log2FC, row.names(DEG_res))
  gsea_res_list[[type]] <- fgsea(pathways = pathway_list, 
                  stats    = ranks,
                  minSize  = 20,
                  maxSize  = 5000) %>% mutate(cellType = type)
}
saveRDS(gsea_res_list, "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_nonTLEControl_glial/Wilcox/gsea_res_list.RDS")

#create tabular format for results 
gsea_res_list <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_nonTLEControl_glial/Wilcox/gsea_res_list.RDS")
mat_for_plot_full <- do.call(rbind, gsea_res_list)
mat_for_plot_full$cellType <- "Glial_Pseudobulk"
mat_for_plot_list[["Glial"]] <- mat_for_plot_full

#####################
#CASE/CONTROL GSEA (NEURON)
#####################
#write out version of DEG files with only protein coding genes
protein_coding_genes <- setNames(unlist(data.table::fread("/n/groups/walsh/indData/Maya/sex_differences_in_the_heart/human_protein_coding_genes.txt", header = TRUE)), NULL)
dir.create("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_NonTLEControl_Neurons/wilcox_10/Protein-Coding/")
for (file in list.files("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_NonTLEControl_Neurons/wilcox_10/All/", pattern = ".RDS"))
{
  res <- readRDS(paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_NonTLEControl_Neurons/wilcox_10/All/", file))
  res <- res %>% filter(row.names(res) %in% protein_coding_genes)
  saveRDS(res, paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_NonTLEControl_Neurons/wilcox_10/Protein-Coding/", file))
  write.csv(res, paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_NonTLEControl_Neurons/wilcox_10/Protein-Coding/", gsub(".RDS", "", file), ".csv"))
}

#get hallmark gene sets
h_gene_sets = msigdbr(species = "human", category = "H")
hallmark_gs <- lapply(unique(h_gene_sets$gs_name), function(x) h_gene_sets$gene_symbol[which(h_gene_sets$gs_name == x)])
names(hallmark_gs) <- unique(h_gene_sets$gs_name)
pathway_list <- hallmark_gs

#add epilepsy gene sets
##DisGeNET
epilepsy_gene_set <- jsonlite::fromJSON(txt="/n/groups/walsh/indData/Maya/FCD_project/raw_data/Ontologies/disgenet_epilepsy_json_string.json")
epilepsy_gene_set <-  epilepsy_gene_set$associations$gene$symbol
pathway_list[["epilepsy_DisGeNET"]] <- epilepsy_gene_set

##Macnee et al., 
macnee_data <- read.csv("/n/groups/walsh/indData/Maya/FCD_project/raw_data/Ontologies/macnee_epilepsy_genes.csv", header = TRUE)
pathway_list[["epilepsy_macnee_high_conf"]] <- (macnee_data %>% filter(Classification == "tier 1"))$Symbol
pathway_list[["epilepsy_macnee"]] <- macnee_data$Symbol

#add synapse gene sets 
##SynGO
syngo_list <- read_excel("/n/groups/walsh/indData/Maya/FCD_project/raw_data/Ontologies/SynGo/syngo_annotations.xlsx") %>% as.data.frame()
syngo_list <- split(syngo_list, syngo_list$go_name)
syngo_list <- lapply(syngo_list, function(x) x$hgnc_symbol)
names(syngo_list) <- sapply(names(syngo_list), function(x) strsplit(x, "(", fixed = TRUE)[[1]][1])
pathway_list <- c(pathway_list , syngo_list) 

#add autism gene set 
sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv", header = TRUE)
SFARI_genes <- unique(sfari$gene.symbol)
high_conf_SFARI_genes <- unique((sfari %>% filter(gene.score == 1))$gene.symbol)
pathway_list[["sfari"]] <- SFARI_genes
pathway_list[["high_conf_sfari"]] <- high_conf_SFARI_genes

#run gsea
wd <- "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_NonTLEControl_Neurons/wilcox_10/Protein-Coding"
gsea_res_list <- list()
for (file in list.files(wd, ".RDS"))
{
  type <- gsub("_DEG_res.RDS", "", file)
  print(type)
  DEG_res <- readRDS(paste0(wd, "/", file)) %>% arrange(desc(avg_log2FC)) %>% filter(pct.1 > 0.10) %>% filter(pct.2 > 0.10)
  ranks <- setNames(DEG_res$avg_log2FC, row.names(DEG_res))
  gsea_res_list[[type]] <- fgsea(pathways = pathway_list, 
                  stats    = ranks,
                  minSize  = 20,
                  maxSize  = 5000) %>% mutate(cellType = type)
}
saveRDS(gsea_res_list, "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_nonTLEControl_neuron/Wilcox/gsea_res_list.RDS")

#create tabular format for results 
gsea_res_list <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_nonTLEControl_neuron/Wilcox/gsea_res_list.RDS")
mat_for_plot_full <- do.call(rbind, gsea_res_list)
mat_for_plot_full$cellType <- "Neuron_Pseudobulk"
mat_for_plot_list[["Neuron"]] <- mat_for_plot_full

#####################
#CASE/CONTROL GSEA (PSEUDOBULK FULL)
#####################
#write out version of DEG files with only protein coding genes
protein_coding_genes <- setNames(unlist(data.table::fread("/n/groups/walsh/indData/Maya/sex_differences_in_the_heart/human_protein_coding_genes.txt", header = TRUE)), NULL)
dir.create("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Pseudobulk/wilcox_10/Protein-Coding/")
for (file in list.files("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Pseudobulk/wilcox_10/All/", pattern = ".RDS"))
{
  res <- readRDS(paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Pseudobulk/wilcox_10/All/", file))
  res <- res %>% filter(row.names(res) %in% protein_coding_genes)
  saveRDS(res, paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Pseudobulk/wilcox_10/Protein-Coding/", file))
  write.csv(res, paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Pseudobulk/wilcox_10/Protein-Coding/", gsub(".RDS", "", file), ".csv"))
}

#get hallmark gene sets
h_gene_sets = msigdbr(species = "human", category = "H")
hallmark_gs <- lapply(unique(h_gene_sets$gs_name), function(x) h_gene_sets$gene_symbol[which(h_gene_sets$gs_name == x)])
names(hallmark_gs) <- unique(h_gene_sets$gs_name)
pathway_list <- hallmark_gs

#add epilepsy gene sets
##DisGeNET
epilepsy_gene_set <- jsonlite::fromJSON(txt="/n/groups/walsh/indData/Maya/FCD_project/raw_data/Ontologies/disgenet_epilepsy_json_string.json")
epilepsy_gene_set <-  epilepsy_gene_set$associations$gene$symbol
pathway_list[["epilepsy_DisGeNET"]] <- epilepsy_gene_set

##Macnee et al., 
macnee_data <- read.csv("/n/groups/walsh/indData/Maya/FCD_project/raw_data/Ontologies/macnee_epilepsy_genes.csv", header = TRUE)
pathway_list[["epilepsy_macnee_high_conf"]] <- (macnee_data %>% filter(Classification == "tier 1"))$Symbol
pathway_list[["epilepsy_macnee"]] <- macnee_data$Symbol

#add synapse gene sets 
##SynGO
syngo_list <- read_excel("/n/groups/walsh/indData/Maya/FCD_project/raw_data/Ontologies/SynGo/syngo_annotations.xlsx") %>% as.data.frame()
syngo_list <- split(syngo_list, syngo_list$go_name)
syngo_list <- lapply(syngo_list, function(x) x$hgnc_symbol)
names(syngo_list) <- sapply(names(syngo_list), function(x) strsplit(x, "(", fixed = TRUE)[[1]][1])
pathway_list <- c(pathway_list , syngo_list) 

#add autism gene set 
sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv", header = TRUE)
SFARI_genes <- unique(sfari$gene.symbol)
high_conf_SFARI_genes <- unique((sfari %>% filter(gene.score == 1))$gene.symbol)
pathway_list[["sfari"]] <- SFARI_genes
pathway_list[["high_conf_sfari"]] <- high_conf_SFARI_genes

#run gsea
wd <- "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Pseudobulk/wilcox_10/Protein-Coding"
gsea_res_list <- list()
for (file in list.files(wd, ".RDS"))
{
  type <- gsub("_DEG_res.RDS", "", file)
  print(type)
  DEG_res <- readRDS(paste0(wd, "/", file)) %>% arrange(desc(avg_log2FC)) %>% filter(pct.1 > 0.10) %>% filter(pct.2 > 0.10)
  ranks <- setNames(DEG_res$avg_log2FC, row.names(DEG_res))
  gsea_res_list[[type]] <- fgsea(pathways = pathway_list, 
                  stats    = ranks,
                  minSize  = 20,
                  maxSize  = 5000) %>% mutate(cellType = type)
}
saveRDS(gsea_res_list, "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_nonTLEControl_pseudobulk/Wilcox/gsea_res_list.RDS")

#create tabular format for results 
gsea_res_list <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_nonTLEControl_pseudobulk/Wilcox/gsea_res_list.RDS")
mat_for_plot_full <- do.call(rbind, gsea_res_list)
mat_for_plot_full$cellType <- "Pseudobulk"
mat_for_plot_list[["Pseudobulk"]] <- mat_for_plot_full

#####################
#PLOT
#####################
rm(mat_for_plot_full)
mat_for_plot_full <- do.call(rbind, mat_for_plot_list)

#add additional formatting changes to tidy things up 
##filter out any pathways that are not significant across all comparisons 
mat_for_plot_full <- mat_for_plot_full %>%
  group_by(pathway) %>%
  filter(any(padj <= 0.05)) %>%
  ungroup()

##convert non significant p-values to a NES of 0  & restrict to only hallmark pathways
mat_for_plot <- mat_for_plot_full %>%
  filter(grepl("HALLMARK", pathway)) %>% as.data.frame() %>%
  mutate(NES = ifelse(padj > 0.05, 0, NES)) %>%
  dplyr::select(cellType, pathway, NES) %>%
  pivot_wider(names_from = cellType, values_from = NES) %>% 
  replace(is.na(.), 0) %>%  # Replace NA values with 0
  column_to_rownames(var = "pathway") %>% 
  as.data.frame()

#create plot
paletteLength <- 50
myBreaks <- c(seq(min(mat_for_plot), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(mat_for_plot)/paletteLength, max(mat_for_plot), length.out=floor(paletteLength/2)))
pheatmap(mat_for_plot, breaks = myBreaks, fontsize_number = 24, cluster_cols = TRUE, color=colorRampPalette(c("navy", "white", "red"))(50),  main = "Hallmark Enriched Pathways", filename = "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_nonTLEControl_pseudobulk/Wilcox/caseVscontrolGSEA_Hallmark.pdf", height = 7, width = 8)

#note: uses matrix generated for panel A
mat_for_plot <- mat_for_plot_full %>%
  filter(!grepl("HALLMARK", pathway)) %>% as.data.frame() %>%
  mutate(NES = ifelse(padj > 0.05, 0, NES)) %>%
  dplyr::select(cellType, pathway, NES) %>%
  pivot_wider(names_from = cellType, values_from = NES) %>% 
  replace(is.na(.), 0) %>%  # Replace NA values with 0
  column_to_rownames(var = "pathway") %>% 
  as.data.frame()

#create plot
paletteLength <- 50
myBreaks <- seq(min(mat_for_plot), 0, length.out = paletteLength + 1)
pheatmap(mat_for_plot, breaks = myBreaks, fontsize_number = 24, cluster_cols = TRUE,  color = colorRampPalette(c("navy", "white"))(paletteLength),  main = "Epilepsy Enriched Pathways", filename = "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_nonTLEControl_pseudobulk/Wilcox/caseVscontrolGSEA_Epilepsy.pdf", height = 6, width = 10)

#####################
#SUPPLEMENTARY TABLE
#####################
DEG_res_list <- list()
DEG_res_list[["Glial"]] <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_NonTLEControl_Glial/wilcox_10/Protein-Coding/DEG_res.RDS")
DEG_res_list[["Neuron"]] <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_NonTLEControl_Neurons/wilcox_10/Protein-Coding/DEG_res.RDS")
DEG_res_list[["Pseudobulk"]] <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Pseudobulk/wilcox_10/Protein-Coding/DEG_res.RDS")
gsea_res_list <- list()
gsea_res_list[["Glial"]] <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_nonTLEControl_glial/Wilcox/gsea_res_list.RDS") %>% as.data.frame() %>% dplyr::select(-DEG_res.RDS.cellType) %>% mutate(cellType = "Glial")
gsea_res_list[["Neuron"]] <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_nonTLEControl_neuron/Wilcox/gsea_res_list.RDS") %>% as.data.frame() %>% dplyr::select(-DEG_res.RDS.cellType) %>% mutate(cellType = "Neuron")
gsea_res_list[["Pseudobulk"]] <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_nonTLEControl_pseudobulk/Wilcox/gsea_res_list.RDS") %>% as.data.frame() %>% dplyr::select(-DEG_res.RDS.cellType) %>% mutate(cellType = "Pseudobulk")
gsea_res_df <- do.call(rbind, gsea_res_list)
colnames(gsea_res_df) <- sapply(colnames(gsea_res_df), function(x) gsub("DEG_res.RDS.", "",x))
row.names(gsea_res_df) <- NULL
DEG_res_list[["GSEA_Results"]] <- gsea_res_df

write.xlsx(DEG_res_list, file = "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_nonTLEControl_pseudobulk/Wilcox/deg_results_case_control_supp_table.xlsx", row.names = TRUE, col.names = TRUE)

