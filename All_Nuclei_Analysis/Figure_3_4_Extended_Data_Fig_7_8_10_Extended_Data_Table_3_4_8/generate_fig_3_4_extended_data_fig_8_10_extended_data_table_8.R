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
library(xlsx)
library(ggrepel)
library(gridExtra)

#source("/n/groups/walsh/indData/Maya/microCHIP_AD_project/GOT/FCD_Analysis/Final_Analysis_For_First_Submission/util.R")

# seurat_obj <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Seurat_Objects/8_processed_seurat_obj_integrated_azimuthAnnotations_ourAnnotations_doubletsFiltered_FCD1Filtered.RDS")
# metadata <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Seurat_Objects/8_processed_seurat_obj_integrated_azimuthAnnotations_ourAnnotations_doubletsFiltered_FCD1Filtered_metadata.RDS")

# #####################
# #CASE/CONTROL GSEA (WILCOX)
# #####################
# #write out version of DEG files with only protein coding genes
# protein_coding_genes <- setNames(unlist(data.table::fread("/n/groups/walsh/indData/Maya/sex_differences_in_the_heart/human_protein_coding_genes.txt", header = TRUE)), NULL)
# dir.create("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_Control/wilcox_10/Protein-Coding/")
# for (file in list.files("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_Control/wilcox_10/All/", pattern = ".RDS"))
# {
# 	res <- readRDS(paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_Control/wilcox_10/All/", file))
# 	res <- res %>% filter(row.names(res) %in% protein_coding_genes)
# 	saveRDS(res, paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_Control/wilcox_10/Protein-Coding/", file))
#   write.csv(res, paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_Control/wilcox_10/Protein-Coding/", gsub(".RDS", "", file), ".csv"))
# }

# #get hallmark gene sets
# h_gene_sets = msigdbr(species = "human", category = "H")
# hallmark_gs <- lapply(unique(h_gene_sets$gs_name), function(x) h_gene_sets$gene_symbol[which(h_gene_sets$gs_name == x)])
# names(hallmark_gs) <- unique(h_gene_sets$gs_name)
# pathway_list <- hallmark_gs

# #add epilepsy gene sets
# ##DisGeNET
# epilepsy_gene_set <- jsonlite::fromJSON(txt="/n/groups/walsh/indData/Maya/FCD_project/raw_data/Ontologies/disgenet_epilepsy_json_string.json")
# epilepsy_gene_set <-  epilepsy_gene_set$associations$gene$symbol
# pathway_list[["epilepsy_DisGeNET"]] <- epilepsy_gene_set

# ##Macnee et al., 
# macnee_data <- read.csv("/n/groups/walsh/indData/Maya/FCD_project/raw_data/Ontologies/macnee_epilepsy_genes.csv", header = TRUE)
# pathway_list[["epilepsy_macnee_high_conf"]] <- (macnee_data %>% filter(Classification == "tier 1"))$Symbol
# pathway_list[["epilepsy_macnee"]] <- macnee_data$Symbol

# #add synapse gene sets 
# ##SynGO
# syngo_list <- read_excel("/n/groups/walsh/indData/Maya/FCD_project/raw_data/Ontologies/SynGo/syngo_annotations.xlsx") %>% as.data.frame()
# syngo_list <- split(syngo_list, syngo_list$go_name)
# syngo_list <- lapply(syngo_list, function(x) x$hgnc_symbol)
# names(syngo_list) <- sapply(names(syngo_list), function(x) strsplit(x, "(", fixed = TRUE)[[1]][1])
# pathway_list <- c(pathway_list , syngo_list) 

# #add autism gene set 
# sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv", header = TRUE)
# SFARI_genes <- unique(sfari$gene.symbol)
# high_conf_SFARI_genes <- unique((sfari %>% filter(gene.score == 1))$gene.symbol)
# pathway_list[["sfari"]] <- SFARI_genes
# pathway_list[["high_conf_sfari"]] <- high_conf_SFARI_genes

# #run gsea
# wd <- "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_Control/wilcox_10/Protein-Coding"
# gsea_res_list <- list()
# for (file in list.files(wd, ".RDS"))
# {
#   type <- gsub("_DEG_res.RDS", "", file)
#   print(type)
#   DEG_res <- readRDS(paste0(wd, "/", file)) %>% arrange(desc(avg_log2FC)) %>% filter(pct.1 > 0.10) %>% filter(pct.2 > 0.10)
#   ranks <- setNames(DEG_res$avg_log2FC, row.names(DEG_res))
#   gsea_res_list[[type]] <- fgsea(pathways = pathway_list, 
#                   stats    = ranks,
#                   minSize  = 20,
#                   maxSize  = 5000) %>% mutate(cellType = type)
# }
# saveRDS(gsea_res_list, "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_control/Wilcox/gsea_res_list.RDS")

# #create tabular format for results 
# gsea_res_list <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_control/Wilcox/gsea_res_list.RDS")
# mat_for_plot_full <- do.call(rbind, gsea_res_list)

# #add additional formatting changes to tidy things up 
# ##filter out any pathways that are not significant across all comparisons 
# mat_for_plot_full <- mat_for_plot_full %>%
#   group_by(pathway) %>%
#   filter(any(padj <= 0.05)) %>%
#   ungroup()

# ##convert non significant p-values to a NES of 0  & restrict to only hallmark pathways
# mat_for_plot <- mat_for_plot_full %>%
#   filter(grepl("HALLMARK", pathway)) %>% as.data.frame() %>%
#   mutate(NES = ifelse(padj > 0.05, 0, NES)) %>%
#   dplyr::select(cellType, pathway, NES) %>%
#   pivot_wider(names_from = cellType, values_from = NES) %>% 
#   replace(is.na(.), 0) %>%  # Replace NA values with 0
#   column_to_rownames(var = "pathway") %>% 
#   as.data.frame()

# #create plot
# paletteLength <- 50
# myBreaks <- c(seq(min(mat_for_plot), 0, length.out=ceiling(paletteLength/2) + 1), 
#               seq(max(mat_for_plot)/paletteLength, max(mat_for_plot), length.out=floor(paletteLength/2)))
# pheatmap(mat_for_plot, breaks = myBreaks, fontsize_number = 24, cluster_cols = TRUE, color=colorRampPalette(c("navy", "white", "red"))(50),  main = "Hallmark Enriched Pathways", filename = "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_control/Wilcox/caseVscontrolGSEA_Hallmark.pdf", height = 7, width = 8)

# #note: uses matrix generated for panel A
# mat_for_plot <- mat_for_plot_full %>%
#   filter(!grepl("HALLMARK", pathway)) %>% as.data.frame() %>%
#   mutate(NES = ifelse(padj > 0.05, 0, NES)) %>%
#   dplyr::select(cellType, pathway, NES) %>%
#   pivot_wider(names_from = cellType, values_from = NES) %>% 
#   replace(is.na(.), 0) %>%  # Replace NA values with 0
#   column_to_rownames(var = "pathway") %>% 
#   as.data.frame()

# #create plot
# paletteLength <- 50
# myBreaks <- c(seq(min(mat_for_plot), 0, length.out=ceiling(paletteLength/2) + 1), 
#               seq(max(mat_for_plot)/paletteLength, max(mat_for_plot), length.out=floor(paletteLength/2)))
# pheatmap(mat_for_plot, breaks = myBreaks, fontsize_number = 24, cluster_cols = TRUE, color=colorRampPalette(c("navy", "white", "red"))(50),  main = "Epilepsy Enriched Pathways", filename = "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_control/Wilcox/caseVscontrolGSEA_Epilepsy.pdf", height = 6, width = 10)

# ######################
# #CASE/CONTROL GSEA (MAST)
# ######################
# #write out version of DEG files with only protein coding genes
# protein_coding_genes <- setNames(unlist(data.table::fread("/n/groups/walsh/indData/Maya/sex_differences_in_the_heart/human_protein_coding_genes.txt", header = TRUE)), NULL)
# dir.create("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_Control/mast_10/Protein-Coding/")
# for (file in list.files("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_Control/mast_10/All/", ".RDS"))
# {
#   res <- readRDS(paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_Control/mast_10/All/", file))
#   res <- res %>% filter(row.names(res) %in% protein_coding_genes)
#   saveRDS(res, paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_Control/mast_10/Protein-Coding/", file))
#   write.csv(res, paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_Control/mast_10/Protein-Coding/", gsub(".RDS", "", file), ".csv"))
# }

# #get hallmark gene sets
# h_gene_sets = msigdbr(species = "human", category = "H")
# hallmark_gs <- lapply(unique(h_gene_sets$gs_name), function(x) h_gene_sets$gene_symbol[which(h_gene_sets$gs_name == x)])
# names(hallmark_gs) <- unique(h_gene_sets$gs_name)
# pathway_list <- hallmark_gs

# #add epilepsy gene sets
# ##DisGeNET
# epilepsy_gene_set <- jsonlite::fromJSON(txt="/n/groups/walsh/indData/Maya/FCD_project/raw_data/Ontologies/disgenet_epilepsy_json_string.json")
# epilepsy_gene_set <-  epilepsy_gene_set$associations$gene$symbol
# pathway_list[["epilepsy_DisGeNET"]] <- epilepsy_gene_set

# ##Macnee et al., 
# macnee_data <- read.csv("/n/groups/walsh/indData/Maya/FCD_project/raw_data/Ontologies/macnee_epilepsy_genes.csv", header = TRUE)
# pathway_list[["epilepsy_macnee_high_conf"]] <- (macnee_data %>% filter(Classification == "tier 1"))$Symbol
# pathway_list[["epilepsy_macnee"]] <- macnee_data$Symbol

# #add synapse gene sets 
# ##SynGO
# syngo_list <- read_excel("/n/groups/walsh/indData/Maya/FCD_project/raw_data/Ontologies/SynGo/syngo_annotations.xlsx") %>% as.data.frame()
# syngo_list <- split(syngo_list, syngo_list$go_name)
# syngo_list <- lapply(syngo_list, function(x) x$hgnc_symbol)
# names(syngo_list) <- sapply(names(syngo_list), function(x) strsplit(x, "(", fixed = TRUE)[[1]][1])
# pathway_list <- c(pathway_list , syngo_list) 

# #add autism gene set 
# sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv", header = TRUE)
# SFARI_genes <- unique(sfari$gene.symbol)
# high_conf_SFARI_genes <- unique((sfari %>% filter(gene.score == 1))$gene.symbol)
# pathway_list[["sfari"]] <- SFARI_genes
# pathway_list[["high_conf_sfari"]] <- high_conf_SFARI_genes

# #run gsea
# wd <- "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_Control/mast_10/Protein-Coding"
# gsea_res_list <- list()
# for (file in list.files(wd, ".RDS"))
# {
#   type <- gsub("_DEG_res.RDS", "", file)
#   print(type)
#   DEG_res <- readRDS(paste0(wd, "/", file)) %>% arrange(desc(avg_log2FC)) %>% filter(pct.1 > 0.10) %>% filter(pct.2 > 0.10)
#   ranks <- setNames(DEG_res$avg_log2FC, row.names(DEG_res))
#   gsea_res_list[[type]] <- fgsea(pathways = pathway_list, 
#                   stats    = ranks,
#                   minSize  = 20,
#                   maxSize  = 5000) %>% mutate(cellType = type)
# }
# saveRDS(gsea_res_list, "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_control/MAST/gsea_res_list.RDS")

# #create tabular format for results 
# gsea_res_list <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_control/MAST/gsea_res_list.RDS")
# mat_for_plot_full <- do.call(rbind, gsea_res_list)

# #add additional formatting changes to tidy things up 
# ##filter out any pathways that are not significant across all comparisons 
# mat_for_plot_full <- mat_for_plot_full %>%
#   group_by(pathway) %>%
#   filter(any(padj <= 0.05)) %>%
#   ungroup()

# ##convert non significant p-values to a NES of 0  & restrict to only hallmark pathways
# mat_for_plot <- mat_for_plot_full %>%
#   filter(grepl("HALLMARK", pathway)) %>% as.data.frame() %>%
#   mutate(NES = ifelse(padj > 0.05, 0, NES)) %>%
#   dplyr::select(cellType, pathway, NES) %>%
#   pivot_wider(names_from = cellType, values_from = NES) %>% 
#   replace(is.na(.), 0) %>%  # Replace NA values with 0
#   column_to_rownames(var = "pathway") %>% 
#   as.data.frame()

# #create plot
# paletteLength <- 50
# myBreaks <- c(seq(min(mat_for_plot), 0, length.out=ceiling(paletteLength/2) + 1), 
#               seq(max(mat_for_plot)/paletteLength, max(mat_for_plot), length.out=floor(paletteLength/2)))
# pheatmap(mat_for_plot, breaks = myBreaks, fontsize_number = 24, cluster_cols = TRUE, color=colorRampPalette(c("navy", "white", "red"))(50),  main = "Hallmark Enriched Pathways", filename = "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_control/MAST/caseVscontrolGSEA_Hallmark.pdf", height = 7, width = 8)

# #note: uses matrix generated for panel A
# mat_for_plot <- mat_for_plot_full %>%
#   filter(!grepl("HALLMARK", pathway)) %>% as.data.frame() %>%
#   mutate(NES = ifelse(padj > 0.05, 0, NES)) %>%
#   dplyr::select(cellType, pathway, NES) %>%
#   pivot_wider(names_from = cellType, values_from = NES) %>% 
#   replace(is.na(.), 0) %>%  # Replace NA values with 0
#   column_to_rownames(var = "pathway") %>% 
#   as.data.frame()

# #create plot
# paletteLength <- 50
# myBreaks <- c(seq(min(mat_for_plot), 0, length.out=ceiling(paletteLength/2) + 1), 
#               seq(max(mat_for_plot)/paletteLength, max(mat_for_plot), length.out=floor(paletteLength/2)))
# pheatmap(mat_for_plot, breaks = myBreaks, fontsize_number = 24, cluster_cols = TRUE, color=colorRampPalette(c("navy", "white", "red"))(50),  main = "Epilepsy Enriched Pathways", filename = "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_control/MAST/caseVscontrolGSEA_Epilepsy.pdf", height = 6, width = 10)

# #######################
# #MICROGLIA STATE ENRICHMENT - WILCOX
# #######################
# #run enrichment analysis
# DEG_res <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_Control/wilcox_10/Protein-Coding/Micro-PVM_DEG_res.RDS") %>% filter(pct.1 > 0.10) %>% filter(pct.2 > 0.10)
# all_marker_genes <- readRDS("/n/groups/walsh/indData/Maya/microCHIP_AD_project/scType/FinalMaterialForPaper/FinalFigures/PanelD/all_marker_genes.RDS")
# all_marker_genes <- all_marker_genes[grepl("Friedman", names(all_marker_genes))]
# universe <- row.names(DEG_res)

# generate_overlap_p_val <- function(list1, list2, universe)
# {
# 	list1 <- intersect(list1, universe)
# 	list2 <- intersect(list2, universe)

# 	inList1AndList2 <- length(intersect(list1, list2))
# 	inList1AndNotList2 <- length(setdiff(list1, list2))
# 	inList2AndNotList1 <- length(setdiff(list2, list1))
# 	inNeither <- length(setdiff(universe, c(list1, list2)))

# 	mat <- matrix(c(inList1AndList2, inList1AndNotList2, inList2AndNotList1, inNeither), nrow = 2)
# 	pval <- fisher.test(mat, alternative = "greater")$p.value

# 	return(pval)
# }
# cluster_marker_enrichment_res <- sapply(all_marker_genes, function(x) generate_overlap_p_val(row.names(DEG_res %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 0.5)), x, universe)) %>% as.data.frame()
# colnames(cluster_marker_enrichment_res) <- "pval"
# cluster_marker_enrichment_res$padj <- p.adjust(cluster_marker_enrichment_res$pval, method = "BH")
# cluster_marker_enrichment_res$PConvert <- -1 * log10(cluster_marker_enrichment_res$padj)
# cluster_marker_enrichment_res$Description <- row.names(cluster_marker_enrichment_res)
# cluster_marker_enrichment_res$num_overlaps <- sapply(cluster_marker_enrichment_res$Description, function(x) length(intersect(row.names(DEG_res %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 0.5)), all_marker_genes[[x]])))
# cluster_marker_enrichment_res <- cluster_marker_enrichment_res %>% arrange(desc(PConvert))
# cluster_marker_enrichment_res$Description <- factor(cluster_marker_enrichment_res$Description, levels = unique(cluster_marker_enrichment_res$Description))
# cluster_marker_enrichment_res$"Significance" <- sapply(cluster_marker_enrichment_res$padj, function(x) ifelse(x > 0.05, "Not Significant", "Significant"))

# #create enrichment plot
# pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/microglia_inflammation_case_vs_control/wilcox/microgliaStateEnrichment.pdf", width = 4, height = 4)
# p <- ggplot(cluster_marker_enrichment_res,aes(x=forcats::fct_rev(Description), y=PConvert,size = num_overlaps, col = Significance))
# p + geom_point() + geom_hline(yintercept=-log10(0.05),linetype="dashed")  + xlab("Microglia State") + ylab("-log10(P-value)")  + theme_classic(base_size = 8) +theme(axis.text=element_text(size=8),axis.title=element_text(size=8)) +  scale_size_continuous(name="Number of DEGs in Microglia State \nMarker Gene List") + scale_color_manual(values=c("blue", "brown1")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()

# #create volcano plot
# vol_plot_list <- list()
# for (state in row.names(cluster_marker_enrichment_res %>% filter(padj < 0.05)))
# #for (state in row.names(cluster_marker_enrichment_res))
# {
#   #get genes to highlight
#   current_marker_genes <- all_marker_genes[[state]]
#   sig_genes <- row.names(DEG_res %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 0.5))
#   sig_genes_overlap_marker_genes_list <- lapply(current_marker_genes, function(x) intersect(x, sig_genes))
#   sig_genes_overlap_marker_genes_list <- unique(do.call(c, sig_genes_overlap_marker_genes_list[order(sapply(sig_genes_overlap_marker_genes_list, length), decreasing = FALSE)]))
#   genes_for_volcano_plot <- sig_genes_overlap_marker_genes_list 

#   #fix p values of 0 
#   DEG_res_og <- DEG_res
#   DEG_res <- DEG_res %>% filter(p_val !=0)

#   DEG_res$symbol <- row.names(DEG_res)
#   DEG_res$Is_Sig_Gene <- sapply(row.names(DEG_res), function(x) ifelse(x %in% genes_for_volcano_plot, "Yes", "No"))
#   DEG_res_is_sig <- DEG_res %>% filter(Is_Sig_Gene == "Yes") 
#   cols <- c("Yes" = "indianred", "No" = "lightgray") 
#   alphas <- c("Yes" = 1, "No" =  0.3)
#   vol_plot <- DEG_res %>%
#   ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj))) + 
#   geom_point(shape = 21, alpha = 0.8, color = "gray") + 
#   geom_point(data = DEG_res_is_sig,
#               shape = 21,
#               fill = "indianred",
#               colour = "black", 
#               alpha = 1)  +   
#   geom_label_repel(data = DEG_res_is_sig, # Add labels last to appear as the top layer  
#                     aes(label = symbol),
#                     size = 1.8,
#                     nudge_y = 1,
#                     min.segment.length = unit(0, 'lines'),
#                     max.overlaps = Inf) +
#   ylab("-Log Adjusted P-Value") + 
#   xlab("Log-2 Fold-Change") + 
#   theme_minimal() + 
#   theme(text=element_text(size=6.5))  + ggtitle(state) +
#     geom_hline(yintercept = -1 *log10(0.05), color = "red", linetype = "dashed") + 
#     geom_vline(xintercept = -0.5, color = "red", linetype = "dashed") + 
#     geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") 

#   #add volcano plot to list
#   vol_plot_list[[state]] <- vol_plot
# }

# pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/microglia_inflammation_case_vs_control/wilcox/microglial_volcano_plots.pdf", height = 4, width = 9)
# print(wrap_plots(vol_plot_list, nrow = 1))
# dev.off()

# #######################
# #MICROGLIA STATE ENRICHMENT - MAST
# #######################
# #run enrichment analysis
# DEG_res <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_Control/mast_10/Protein-Coding/Micro-PVM_DEG_res.RDS") %>% filter(pct.1 > 0.10) %>% filter(pct.2 > 0.10)
# all_marker_genes <- readRDS("/n/groups/walsh/indData/Maya/microCHIP_AD_project/scType/FinalMaterialForPaper/FinalFigures/PanelD/all_marker_genes.RDS")
# all_marker_genes <- all_marker_genes[grepl("Friedman", names(all_marker_genes))]
# universe <- row.names(DEG_res)

# generate_overlap_p_val <- function(list1, list2, universe)
# {
# 	list1 <- intersect(list1, universe)
# 	list2 <- intersect(list2, universe)

# 	inList1AndList2 <- length(intersect(list1, list2))
# 	inList1AndNotList2 <- length(setdiff(list1, list2))
# 	inList2AndNotList1 <- length(setdiff(list2, list1))
# 	inNeither <- length(setdiff(universe, c(list1, list2)))

# 	mat <- matrix(c(inList1AndList2, inList1AndNotList2, inList2AndNotList1, inNeither), nrow = 2)
# 	pval <- fisher.test(mat, alternative = "greater")$p.value

# 	return(pval)
# }
# cluster_marker_enrichment_res <- sapply(all_marker_genes, function(x) generate_overlap_p_val(row.names(DEG_res %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 0.5)), x, universe)) %>% as.data.frame()
# colnames(cluster_marker_enrichment_res) <- "pval"
# cluster_marker_enrichment_res$padj <- p.adjust(cluster_marker_enrichment_res$pval, method = "BH")
# cluster_marker_enrichment_res$PConvert <- -1 * log10(cluster_marker_enrichment_res$padj)
# cluster_marker_enrichment_res$Description <- row.names(cluster_marker_enrichment_res)
# cluster_marker_enrichment_res$num_overlaps <- sapply(cluster_marker_enrichment_res$Description, function(x) length(intersect(row.names(DEG_res %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 0.5)), all_marker_genes[[x]])))
# cluster_marker_enrichment_res <- cluster_marker_enrichment_res %>% arrange(desc(PConvert))
# cluster_marker_enrichment_res$Description <- factor(cluster_marker_enrichment_res$Description, levels = unique(cluster_marker_enrichment_res$Description))
# cluster_marker_enrichment_res$"Significance" <- sapply(cluster_marker_enrichment_res$padj, function(x) ifelse(x > 0.05, "Not Significant", "Significant"))

# #create enrichment plot
# pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/microglia_inflammation_case_vs_control/MAST/microgliaStateEnrichment.pdf", width = 4, height = 4)
# p <- ggplot(cluster_marker_enrichment_res,aes(x=forcats::fct_rev(Description), y=PConvert,size = num_overlaps, col = Significance))
# p + geom_point() + geom_hline(yintercept=-log10(0.05),linetype="dashed")  + xlab("Microglia State") + ylab("-log10(P-value)")  + theme_classic(base_size = 8) +theme(axis.text=element_text(size=8),axis.title=element_text(size=8)) +  scale_size_continuous(name="Number of DEGs in Microglia State \nMarker Gene List") + scale_color_manual(values=c("blue", "brown1")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()

# #create volcano plot
# vol_plot_list <- list()
# #for (state in row.names(cluster_marker_enrichment_res %>% filter(padj < 0.05)))
# for (state in row.names(cluster_marker_enrichment_res))
# {
#   #get genes to highlight
#   current_marker_genes <- all_marker_genes[[state]]
#   sig_genes <- row.names(DEG_res %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 0.5))
#   sig_genes_overlap_marker_genes_list <- lapply(current_marker_genes, function(x) intersect(x, sig_genes))
#   sig_genes_overlap_marker_genes_list <- unique(do.call(c, sig_genes_overlap_marker_genes_list[order(sapply(sig_genes_overlap_marker_genes_list, length), decreasing = FALSE)]))
#   genes_for_volcano_plot <- sig_genes_overlap_marker_genes_list 

#   #fix p values of 0 
#   DEG_res_og <- DEG_res
#   DEG_res <- DEG_res %>% filter(p_val !=0)

#   DEG_res$symbol <- row.names(DEG_res)
#   DEG_res$Is_Sig_Gene <- sapply(row.names(DEG_res), function(x) ifelse(x %in% genes_for_volcano_plot, "Yes", "No"))
#   DEG_res_is_sig <- DEG_res %>% filter(Is_Sig_Gene == "Yes") 
#   cols <- c("Yes" = "indianred", "No" = "lightgray") 
#   alphas <- c("Yes" = 1, "No" =  0.3)
#   vol_plot <- DEG_res %>%
#   ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj))) + 
#   geom_point(shape = 21, alpha = 0.8, color = "gray") + 
#   geom_point(data = DEG_res_is_sig,
#               shape = 21,
#               fill = "indianred",
#               colour = "black", 
#               alpha = 1)  +   
#   geom_label_repel(data = DEG_res_is_sig, # Add labels last to appear as the top layer  
#                     aes(label = symbol),
#                     size = 1.8,
#                     nudge_y = 1,
#                     min.segment.length = unit(0, 'lines'),
#                     max.overlaps = Inf) +
#   ylab("-Log Adjusted P-Value") + 
#   xlab("Log-2 Fold-Change") + 
#   theme_minimal() + 
#   theme(text=element_text(size=6.5))  + ggtitle(state) +
#     geom_hline(yintercept = -1 *log10(0.05), color = "red", linetype = "dashed") + 
#     geom_vline(xintercept = -0.5, color = "red", linetype = "dashed") + 
#     geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") 

#   #add volcano plot to list
#   vol_plot_list[[state]] <- vol_plot
# }

# pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/microglia_inflammation_case_vs_control/MAST/microglial_volcano_plots.pdf", height = 4, width = 9)
# print(wrap_plots(vol_plot_list, nrow = 1))
# dev.off()

#######################
#CELL CHAT
#######################
# #add in higher level annots
# higher_level_annots <- c(
# "Oligo" = "Oligo/OPC", 
# "OPC" = "Oligo/OPC",
# "Astro" = "Astro",
# "L2/3 IT" = "UL_ExN",
# "Micro-PVM" = "Microglia", 
# "Pvalb" = "IN_MGE", 
# "L4 IT" = "DL_ExN", 
# "L5 IT" = "DL_ExN", 
# "L5 ET" = "DL_ExN", 
# "Vip" = "IN_CGE", 
# "Sst" = "IN_MGE", 
# "Sncg" = "IN_CGE",
# "Lamp5" = "IN_CGE", 
# "L5/6 NP" = "DL_ExN", 
# "L6 CT" = "DL_ExN", 
# "L6 IT" = "DL_ExN", 
# "L6 IT Car3" = "DL_ExN", 
# "L6b" = "DL_ExN", 
# "Endo" = "Endo", 
# "VLMC" = "VLMC"
# )
# seurat_obj$final_class <- setNames(higher_level_annots[seurat_obj$final_annots], names(seurat_obj$final_annots))

# seurat_obj$samples <- seurat_obj$cleaned_ident
# seurat_obj_og <- seurat_obj
# Idents(seurat_obj_og) <- "status"

# #create cellchat object for cases 
# seurat_obj <- subset(seurat_obj_og, idents = "case")
# cellchat <- createCellChat(object = seurat_obj, group.by = "final_class", assay = "RNA")

# CellChatDB <- CellChatDB.human 
# cellchat@DB <- CellChatDB
# cellchat <- subsetData(cellchat) 

# cellchat <- identifyOverExpressedGenes(cellchat)
# cellchat <- identifyOverExpressedInteractions(cellchat)
# cellchat <- computeCommunProb(cellchat, type = "triMean")
# cellchat <- filterCommunication(cellchat, min.cells = 10)
# cellchat <- computeCommunProbPathway(cellchat)
# cellchat <- aggregateNet(cellchat)
# cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
# saveRDS(cellchat, "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/CellChat_Case_Control/cell_chat_case_collapsed.RDS")

# #create cellchat object for controls 
# seurat_obj <- subset(seurat_obj_og, idents = "control")
# rm(seurat_obj_og)
# cellchat <- createCellChat(object = seurat_obj, group.by = "final_class", assay = "RNA")

# CellChatDB <- CellChatDB.human 
# cellchat@DB <- CellChatDB
# cellchat <- subsetData(cellchat) 

# cellchat <- identifyOverExpressedGenes(cellchat)
# cellchat <- identifyOverExpressedInteractions(cellchat)
# cellchat <- computeCommunProb(cellchat, type = "triMean")
# cellchat <- filterCommunication(cellchat, min.cells = 10)
# cellchat <- computeCommunProbPathway(cellchat)
# cellchat <- aggregateNet(cellchat)
# cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
# saveRDS(cellchat, "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/CellChat_Case_Control/cell_chat_control_collapsed.RDS")
# rm(seurat_obj)
# rm(seurat_obj_og)

#compare cases vs controls
##****************merge the cellchat objects
cellchat_case <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/CellChat_Case_Control/cell_chat_case_collapsed.RDS")
cellchat_control <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/CellChat_Case_Control/cell_chat_control_collapsed.RDS")
object.list <- list(control = cellchat_control, case = cellchat_case)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
all_signaling_pathways <- unique(CellChatDB.human$interaction$pathway)
all_signaling_pathways <- all_signaling_pathways[order(all_signaling_pathways)]

##****************compare number of interactions 
pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/cellchat_case_vs_control/0_cellChat_numInteractions.pdf", width = 12, height = 6)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
print(gg1 + gg2)
dev.off()

##****************create a heatmap showing differential interactions --> in general, excitatory neurons seem to be doing a lot more signaling in cases versus controls 
pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/cellchat_case_vs_control/1_cellChat_differentialHeatmap.pdf", width = 12, height = 6)
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
print(gg1 + gg2)
dev.off()

##****************create a chord plot showing differential iteractions between higher groupings of celltypes
pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/cellchat_case_vs_control/2_cellChat_differentialCircosPlot.pdf", width = 12, height = 6)  
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

#****************identify the signaling changes of specific cell populations
plot_list <- list()
cell_types <- c("DL_ExN", "UL_ExN", "IN_MGE", "IN_CGE", "Microglia", "Oligo/OPC", "Astro", "VLMC", "Endo")
for (type in cell_types) {
    print(type)
    tryCatch({
        plot_list[[type]] <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = type)
    }, error = function(e) {
        message("Error in netAnalysis_signalingChanges_scatter for type: ", type, "\n", e)
    })
}

pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/cellchat_case_vs_control/3_SignalingChangesPerCellType.pdf", width = 16, height = 16)
patchwork::wrap_plots(plot_list, nrow = 3)
dev.off()

##****************identify altered signaling with distinct network architecture and interaction strength
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")

##distance between pathways 
pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/cellchat_case_vs_control/4_overallDistanceBetweenPathways.pdf", width = 16, height = 16)
rankSimilarity(cellchat, type = "functional")
dev.off()

##altered information flow
pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/cellchat_case_vs_control/5_overallChangesInInfoFlow.pdf", width = 8, height = 16)
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
print(gg1)
dev.off()

##****************compare cell-cell communication between cases and controls
pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/cellchat_case_vs_control/6_pathway_activity_btwn_cases_and_controls.pdf", width = 12, height = 6)
all_pathways <- cellchat@netP$control$pathway
for (pathway in all_pathways) {
  tryCatch({
    print(pathway)
    pathways.show <- pathway
    par(mfrow = c(1, 2), xpd = TRUE)
    ht <- list()
    
    for (i in 1:length(object.list)) {
      ht[[i]] <- netVisual_heatmap(
        object.list[[i]], 
        signaling = pathways.show, 
        color.heatmap = "Reds",
        title.name = paste(pathways.show, "signaling", names(object.list)[i])
      )
    }
    
    ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
  }, error = function(e) {
    message("Error encountered for pathway: ", pathway)
    message("Error message: ", e$message)
  })
}
dev.off()

##****************identify dysfunctional signaling by using differential expression analysis 
pos.dataset = "case"
features.name = paste0(pos.dataset, ".merged")
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
net.up <- subsetCommunication(cellchat, net = net, datasets = "case",ligand.logFC = 0.05, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "control",ligand.logFC = -0.05, receptor.logFC = NULL)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
write.table(net.up, "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/cellchat_case_vs_control/7_cellchat_net_up.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(net.down, "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/cellchat_case_vs_control/7_cellchat_net_down.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

##****************compare number of interactions 
ctypes <- c("Astro", "DL_ExN", "Endo", "IN_CGE", "IN_MGE", "Microglia", 
                "Oligo/OPC", "UL_ExN", "VLMC")
pairs_to_test <- list(c("DL_ExN", "IN_CGE"), c("DL_ExN", "IN_MGE"), c("DL_ExN", "UL_ExN"), c("IN_MGE", "IN_MGE"), c("UL_ExN", "UL_ExN"))
pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/cellchat_case_vs_control/8_cellchat_diff_signaling.pdf", width = 8, height = 16)
for (pair in pairs_to_test)
{
  source_idx <- which(ctypes == pair[1])
  target_idx <- which(ctypes == pair[2])

  gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = source_idx, targets.use = target_idx, stacked = T, do.stat = TRUE)
  gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = source_idx, targets.use = target_idx, stacked = F, do.stat = TRUE)
  print(grid.arrange(gg1, gg2, ncol = 2, top = paste("Comparison: ", pair[1], " vs ", pair[2])))
}
dev.off()

