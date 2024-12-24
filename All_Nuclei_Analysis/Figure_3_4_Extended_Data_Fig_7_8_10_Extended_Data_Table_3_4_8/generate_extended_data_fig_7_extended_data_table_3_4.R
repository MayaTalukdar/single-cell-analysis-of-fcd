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

#######################
#SUPPLEMENTARY TABLES - WILCOX
#######################
all_files <- list.files("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_Control/wilcox_10/Protein-Coding/", pattern = ".RDS")
DEG_res_list <- lapply(all_files, function(x) readRDS(paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_Control/wilcox_10/Protein-Coding/", x)) %>% arrange(desc(avg_log2FC)) %>% filter(pct.1 > 0.10) %>% filter(pct.2 > 0.10))
names(DEG_res_list) <- sapply(all_files, function(x) gsub("_DEG_res.RDS", "", x))
DEG_res_list[["GSEA_Results"]] <- do.call(rbind, readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_control/Wilcox/gsea_res_list.RDS")) %>% as.data.frame()
write.xlsx(DEG_res_list, file = "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_control/Wilcox/deg_results_case_control_supp_table.xlsx", row.names = TRUE, col.names = TRUE)

#######################
#SUPPLEMENTARY TABLES - MAST
#######################
all_files <- list.files("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_Control/mast_10/Protein-Coding/", pattern = ".RDS")
DEG_res_list <- lapply(all_files, function(x) readRDS(paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_Control/mast_10/Protein-Coding/", x)) %>% arrange(desc(avg_log2FC)) %>% filter(pct.1 > 0.10) %>% filter(pct.2 > 0.10))
names(DEG_res_list) <- sapply(all_files, function(x) gsub("_DEG_res.RDS", "", x))
DEG_res_list[["GSEA_Results"]] <- do.call(rbind, readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_control/MAST/gsea_res_list.RDS")) %>% as.data.frame()
write.xlsx(DEG_res_list, file = "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_control/MAST/deg_results_case_control_supp_table.xlsx", row.names = TRUE, col.names = TRUE)

#######################
#VOLCANO PLOTS - WILCOX
#######################
vol_plot_list <- list()
for (file in list.files("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_Control/wilcox_10/Protein-Coding/", pattern = ".RDS")) {
  DEG_res <- readRDS(paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/DEG_Results_Case_Control/wilcox_10/Protein-Coding/", file))
  type <- gsub("_DEG_res.RDS", "", file)
  
  # Filter p-values and set thresholds
  DEG_res <- DEG_res %>% filter(p_val_adj != 0) %>% filter(pct.1 > 0.10) %>% filter(pct.2 > 0.10)
  DEG_res$symbol <- row.names(DEG_res)
  
  # Define significance thresholds
  log2FC_threshold <- 0.5
  p_adj_threshold <- 0.05
  
  DEG_res <- DEG_res %>%
    mutate(
      significance = case_when(
        avg_log2FC > log2FC_threshold & p_val_adj < p_adj_threshold ~ "Upregulated",
        avg_log2FC < -log2FC_threshold & p_val_adj < p_adj_threshold ~ "Downregulated",
        TRUE ~ "Not Significant"
      )
    )
  
  # Define colors for significance
  cols <- c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")
  
  # Create volcano plot
  vol_plot <- DEG_res %>%
    ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = cols) +
    ylab("-Log10 Adjusted P-Value") +
    xlab("Log2 Fold-Change") +
    geom_hline(yintercept = -1 *log10(0.05), color = "red", linetype = "dashed") + 
    geom_vline(xintercept = -0.5, color = "red", linetype = "dashed") + 
    geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") + 
    theme(text = element_text(size = 10)) +
    ggtitle(type) + theme_minimal()

  
  # Add volcano plot to list
  vol_plot_list[[type]] <- vol_plot
}

pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/diff_expr_case_vs_control/Wilcox/celltype_volcano_plots.pdf", height = 20, width = 12)
print(wrap_plots(vol_plot_list, ncol = 3))
dev.off()
