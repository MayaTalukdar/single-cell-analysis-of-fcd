library(Seurat)
library(tidyverse)
library(scCustomize)
library(cowplot)
library(patchwork)
library(RColorBrewer)

###########################
#I/O
###########################
# seurat_obj <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Seurat_Objects/8_processed_seurat_obj_integrated_azimuthAnnotations_ourAnnotations_doubletsFiltered_FCD1Filtered.RDS")
# full_metadata <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Seurat_Objects/8_processed_seurat_obj_integrated_azimuthAnnotations_ourAnnotations_doubletsFiltered_FCD1Filtered_metadata.RDS")

# #create dataframes of genotyping results 
# samples <- c("e174", "FC5801", "FC5501")
# metadata_list <- list()
# for (current_sample in samples)
# {
#     print("*********************************")
#     print(paste0("STARTING SAMPLE ", current_sample, "!"))
#     print("*********************************")

#     #read in genotyping results 
#     collapsed_10X_UMI_genotype_results <- read.table(paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Genotyping/", current_sample, "/", "scMosaicHunter_results_", current_sample, "_using_collapsed_10X_umi_with_custom_VAF_strictGenotyping_just_relevant_cols.txt"), header = TRUE)
#     table(collapsed_10X_UMI_genotype_results$genotype)
#     print("*********************************")
#     print(paste0(nrow(collapsed_10X_UMI_genotype_results), " NUCLEI GENOTYPED!"))
#     print("*********************************")

#     #merge with metadata
#     metadata <- full_metadata %>% filter(orig.ident == current_sample)
#     metadata$orig_row_names <- row.names(metadata)
#     metadata$cleaned_cell_name <- sapply(row.names(metadata), function(x) strsplit(x, "-")[[1]][2])
#     print("*********************************")
#     print(paste0(length(which(metadata$cleaned_cell_name %in% collapsed_10X_UMI_genotype_results$best.bc)), " OF THESE NUCLEI WERE INCLUDED IN FINAL SEURAT OBJECT!")) 
#     print("*********************************")

#     orig_n_rows <- nrow(metadata)
#     metadata <- merge(metadata, collapsed_10X_UMI_genotype_results, by.x = "cleaned_cell_name",  by.y = "best.bc", all.x = TRUE)
#     print(orig_n_rows == nrow(metadata))
#     metadata$genotype[which(is.na(metadata$genotype))] <- "NotGenotyped"
#     row.names(metadata) <- metadata$orig_row_names
#     metadata_list[[current_sample]] <- metadata
# }

# #add them back to the seurat object

# metadata_list_merged <- do.call(rbind, metadata_list)
# row.names(metadata_list_merged) <- metadata_list_merged$orig_row_names
# genotype_mapper_vec <- setNames(metadata_list_merged$genotype, row.names(metadata_list_merged))
# table(names(genotype_mapper_vec) %in% row.names(full_metadata))
# seurat_obj$genotype <- setNames(genotype_mapper_vec[row.names(seurat_obj@meta.data)], row.names(seurat_obj@meta.data))
# print(sum(table(seurat_obj$genotype)) == length(genotype_mapper_vec))
# seurat_obj$genotype[which(is.na(seurat_obj$genotype))] <- "CellInNotGenotypedSample"

# #look at breakdown of genotype over a few covariates of interest
# table(seurat_obj$final_annots, seurat_obj$genotype)[,c("hetero", "ref-hom")]
# table(seurat_obj$cleaned_ident, seurat_obj$genotype)[,c("hetero", "ref-hom")]
# saveRDS(seurat_obj, "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Seurat_Objects/9_processed_seurat_obj_integrated_azimuthAnnotations_ourAnnotations_doubletsFiltered_FCD1Filtered_strictGenotyping.RDS")
# saveRDS(seurat_obj@meta.data, "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Seurat_Objects/9_processed_seurat_obj_integrated_azimuthAnnotations_ourAnnotations_doubletsFiltered_FCD1Filtered_strictGenotyping_metadata.RDS")

seurat_obj <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Seurat_Objects/9_processed_seurat_obj_integrated_azimuthAnnotations_ourAnnotations_doubletsFiltered_FCD1Filtered_strictGenotyping.RDS")
Idents(seurat_obj) <- "orig.ident"

#add in final class categories
higher_level_annots <- c(
"Oligo" = "Oligo/OPC", 
"OPC" = "Oligo/OPC",
"Astro" = "Astro",
"L2/3 IT" = "UL_ExN",
"Micro-PVM" = "Microglia", 
"Pvalb" = "IN_MGE", 
"L4 IT" = "DL_ExN", 
"L5 IT" = "DL_ExN", 
"L5 ET" = "DL_ExN", 
"Vip" = "IN_CGE", 
"Sst" = "IN_MGE", 
"Sncg" = "IN_CGE",
"Lamp5" = "IN_CGE", 
"L5/6 NP" = "DL_ExN", 
"L6 CT" = "DL_ExN", 
"L6 IT" = "DL_ExN", 
"L6 IT Car3" = "DL_ExN", 
"L6b" = "DL_ExN", 
"Endo" = "Endo", 
"VLMC" = "VLMC"
)
seurat_obj$final_class <- setNames(higher_level_annots[seurat_obj$final_annots], names(seurat_obj$final_annots))

#add in final supraclass category 
seurat_obj$final_supraclass <- sapply(seurat_obj$final_class, function(x) ifelse(grepl("ExN", x), "Excitatory", ifelse(grepl("IN", x), "Inhibitory", "Non-Neuronal")))

#subset to only genotyped sample
Idents(seurat_obj) <- "orig.ident"
seurat_obj <- subset(seurat_obj, idents = c("e174", "FC5501", "FC5801"))

######################
#UMAP BY CELL TYPE
######################
umap_by_type <- DimPlot(seurat_obj, group.by = "final_annots", reduction = "umap.harmony", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Cell Type")  

######################
#UMAP OF EXPRESSION
######################
pik3ca_seurat_obj <- subset(seurat_obj, idents = "e174")
pik3ca_umap_by_expression <- FeaturePlot_scCustom(pik3ca_seurat_obj, features = "PIK3CA", reduction = "umap.harmony") + NoLegend() + ggtitle("PIK3CA Mut: PIK3CA Expression")
pik3ca_umap_by_expression_with_legend <- FeaturePlot_scCustom(pik3ca_seurat_obj, features = "PIK3CA", reduction = "umap.harmony") + ggtitle("PIK3CA Mut: PIK3CA Expression")

mtor_seurat_obj <- subset(seurat_obj, idents = c("FC5501", "FC5801"))
mtor_umap_by_expression <- FeaturePlot_scCustom(mtor_seurat_obj, features = "MTOR", reduction = "umap.harmony") + NoLegend() + ggtitle("MTOR Mut: MTOR Expression")
mtor_umap_by_expression_with_legend <- FeaturePlot_scCustom(mtor_seurat_obj, features = "MTOR", reduction = "umap.harmony") + ggtitle("MTOR Mut: MTOR Expression")

######################
#UMAP OF GENOTYPING
######################
umap_by_genotyping <- DimPlot(seurat_obj, group.by = "genotype", cols = c("red", "lightgray", "#097969"), reduction = "umap.harmony") + ggtitle("Genotyping")  + NoLegend()
umap_by_genotyping[[1]]$layers[[1]]$aes_params$alpha =  ifelse (seurat_obj@meta.data$genotype == "NotGenotyped", 0.2, 1)

######################
#BARPLOT OF GENOTYPING PER CELL TYPE 
######################
metadata <- seurat_obj@meta.data %>% filter(genotype != "NotGenotyped")
num_cells_df <- table(metadata$genotype, metadata$final_annots) %>% t() %>% as.data.frame()
colnames(num_cells_df) <- c('final_annots', 'genotype', 'numCells')
ordered_levels <- num_cells_df %>%
  group_by(final_annots) %>%
  summarise(total_cells = sum(numCells)) %>%
  arrange(desc(total_cells)) %>%
  pull(final_annots)
num_cells_df$final_annots <- factor(num_cells_df$final_annots, levels = ordered_levels)
barplot <- ggplot(num_cells_df, aes(x = final_annots, y = numCells, fill = genotype)) +
  geom_bar(stat = "identity") +
  labs(x = "final_annots", y = "numCells") +
  scale_fill_manual(values = c("hetero" = "red", "ref-hom" = "#097969")) +  # Updated colors
  theme_minimal()+ ggtitle("Genotyping Per Cell Type") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + NoLegend()

######################
#BARPLOT OF PROPORTION OF MUTANT CELLS - STRICT
######################
metadata <- seurat_obj@meta.data %>% filter(genotype != "NotGenotyped")
pval_list <- list()
prop_list <- list()
for (type in unique((metadata %>% filter(genotype == "hetero"))$final_annots))
{
  #create matrix of num hetero vs ref homo cells
  num_het_cells_type <- nrow(metadata %>% filter(final_annots == type) %>% filter(genotype == "hetero"))
  num_ref_cells_type <- nrow(metadata %>% filter(final_annots == type) %>% filter(genotype == "ref-hom"))
  num_het_cells <- nrow(metadata %>% filter(genotype == "hetero"))
  num_ref_cells <- nrow(metadata %>% filter(genotype == "ref-hom"))
  mat <- matrix(data = c(num_het_cells_type, num_ref_cells_type, num_het_cells, num_ref_cells), nrow = 2)
  pval_list[[type]] <- fisher.test(mat, alternative = "greater")$p.value
  prop_list[[type]] <- mat[1,1]/sum(mat[,1])
}

prop_df <- prop_list %>% as.data.frame() %>% t() %>% as.data.frame()
row.names(prop_df) <- unique((metadata %>% filter(genotype == "hetero"))$final_annots)
colnames(prop_df) <- "prop"
prop_df$type <- row.names(prop_df)
prop_df$type <- factor(prop_df$type, levels = prop_df$type[order(prop_df$prop)])

barplot_mutProp_strict <- ggplot(prop_df, aes(y = type, x = prop, fill = "red")) +
  geom_bar(stat = "identity") +
  labs(y = "type", x = "propMutantCells") +
  theme_minimal()+ ggtitle("Proportion of Mutant Cells") + NoLegend()

write.table(pval_list %>% as.data.frame() %>% t() %>% as.data.frame(), "/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/genotyping_figure_strict/pval_list.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

######################
#BARPLOT OF PROPORTION OF MUTANT CELLS - LOOSE
######################
cell_types_to_keep <- c("Oligo", "Astro", "L6 CT", "L6b", "L2/3 IT", "L4 IT", "L5 IT", "Sst", "Vip", "Sncg")
metadata <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Seurat_Objects/9_processed_seurat_obj_integrated_azimuthAnnotations_ourAnnotations_doubletsFiltered_FCD1Filtered_looseGenotyping_metadata.RDS") %>% filter(genotype != "NotGenotyped")
prop_list <- list()
for (type in unique((metadata %>% filter(genotype == "hetero"))$final_annots))
{
  #create matrix of num hetero vs ref homo cells
  num_het_cells_type <- nrow(metadata %>% filter(final_annots == type) %>% filter(genotype == "hetero"))
  num_ref_cells_type <- nrow(metadata %>% filter(final_annots == type) %>% filter(genotype == "ref-hom"))
  num_het_cells <- nrow(metadata %>% filter(genotype == "hetero"))
  num_ref_cells <- nrow(metadata %>% filter(genotype == "ref-hom"))
  mat <- matrix(data = c(num_het_cells_type, num_ref_cells_type, num_het_cells, num_ref_cells), nrow = 2)
  prop_list[[type]] <- mat[1,1]/sum(mat[,1])
}
prop_list <- prop_list[cell_types_to_keep]
prop_df <- prop_list %>% as.data.frame() %>% t() %>% as.data.frame()
row.names(prop_df) <- cell_types_to_keep
colnames(prop_df) <- "prop"
prop_df$type <- row.names(prop_df)
prop_df$type <- factor(prop_df$type, levels = prop_df$type[order(prop_df$prop)])

barplot_mutProp <- ggplot(prop_df, aes(y = type, x = prop, fill = "red")) +
  geom_bar(stat = "identity") +
  labs(y = "type", x = "propMutantCells") +
  theme_minimal()+ ggtitle("Proportion of Mutant Cells") + NoLegend()

pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/genotyping_figure_strict/barplots_based_on_loose_results.pdf")
print(barplot_mutProp)
dev.off()


######################
#PIE CHART 
######################
pie_df <- seurat_obj@meta.data %>% filter(genotype == "hetero")
pie_df <- table(pie_df$final_supraclass)/sum(table(pie_df$final_supraclass))
pie_df <- as.data.frame(pie_df)
colnames(pie_df) <- c("l1", "prop")
colors <- brewer.pal(n = nrow(pie_df), name = "Set1")
pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/genotyping_figure_strict/pie_chart.pdf", width = 6, height = 6)
colnames(pie_df) <- c("l1", "prop")
ggplot(pie_df, aes(x = "", y = prop, fill = l1)) +
geom_bar(stat = "identity", width = 1) +
coord_polar(theta = "y") +
scale_fill_manual(values = colors) +
theme_void()
dev.off()

######################
#FINAL PLOTS
######################
pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/genotyping_figure_strict/genotyping_figure_strict.pdf", width = 4, height =18)
print(wrap_plots(umap_by_type, pik3ca_umap_by_expression, mtor_umap_by_expression, umap_by_genotyping, barplot, barplot_mutProp_strict, ncol = 1))
dev.off()

######################
#MODELING
######################
#PIK3CA
#set up design matrix and fit model 
gene_to_genotype <- "PIK3CA"
design_mat <- cbind(seurat_obj$final_annots, GetAssayData(seurat_obj)[gene_to_genotype,], seurat_obj$genotype) %>% as.data.frame()
colnames(design_mat) <- c("final_annots", "gene_expr", "genotype")
design_mat$wasGenotyped <- sapply(design_mat$genotype, function(x) ifelse(x != "NotGenotyped", 1, 0))
design_mat$gene_expr <- as.numeric(design_mat$gene_expr)
design_mat$final_annots <- as.factor(design_mat$final_annots)
model <- glm(wasGenotyped ~ final_annots + gene_expr, data = design_mat)

#question #1: what factors affect genotyping probability?
summary(model)

#plot
ready1=data.frame(coef(summary(model)))[-1,]
colnames(ready1)=c("Estimate","SE","Tvalue","Pvalue")
ready1$Coefficient=factor(rownames(ready1),levels=rev(rownames(ready1)))
ready1$LowCI=ready1$Estimate-ready1$SE*2
ready1$HighCI=ready1$Estimate+ready1$SE*2
ready1$Coefficient <- sapply(row.names(ready1), function(x) gsub("final_annots", "", x))
ready1$Coefficient[nrow(ready1)] <- "PIK3CA Expr."
ready1$Coefficient <- as.factor(ready1$Coefficient)
ready1$Coefficient <- relevel(ready1$Coefficient, ref = "PIK3CA Expr.")
ready1$Coefficient <- fct_rev(ready1$Coefficient)
pik3ca_geno_plot <- ggplot(ready1,aes(x=Coefficient,y=Estimate))
pik3ca_geno_plot <- pik3ca_geno_plot + geom_errorbar(aes(ymin=LowCI,ymax=HighCI),width=0.2) + geom_point(shape=16,size=3,fill="white") + geom_hline(yintercept=0,linetype="dashed") + scale_fill_discrete(name="") + xlab("") + ylab("Estimate") + coord_flip() + theme_classic() + theme(text=element_text(size=12)) + ggtitle("PIK3CA")

pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/genotyping_figure_strict/pik3ca_genotyping_eff_modeling.pdf", width = 10, height =4)
print(pik3ca_geno_plot)
dev.off()

#MTOR
#set up design matrix and fit model 
gene_to_genotype <- "MTOR"
design_mat <- cbind(seurat_obj$final_annots, GetAssayData(seurat_obj)[gene_to_genotype,], seurat_obj$genotype) %>% as.data.frame()
colnames(design_mat) <- c("final_annots", "gene_expr", "genotype")
design_mat$wasGenotyped <- sapply(design_mat$genotype, function(x) ifelse(x != "NotGenotyped", 1, 0))
design_mat$gene_expr <- as.numeric(design_mat$gene_expr)
design_mat$final_annots <- as.factor(design_mat$final_annots)
model <- glm(wasGenotyped ~ final_annots + gene_expr, data = design_mat)

#question #1: what factors affect genotyping probability?
summary(model)

#plot
ready1=data.frame(coef(summary(model)))[-1,]
colnames(ready1)=c("Estimate","SE","Tvalue","Pvalue")
ready1$Coefficient=factor(rownames(ready1),levels=rev(rownames(ready1)))
ready1$LowCI=ready1$Estimate-ready1$SE*2
ready1$HighCI=ready1$Estimate+ready1$SE*2
ready1$Coefficient <- sapply(row.names(ready1), function(x) gsub("final_annots", "", x))
ready1$Coefficient[nrow(ready1)] <- "MTOR Expr."
ready1$Coefficient <- as.factor(ready1$Coefficient)
ready1$Coefficient <- relevel(ready1$Coefficient, ref = "MTOR Expr.")
ready1$Coefficient <- fct_rev(ready1$Coefficient)
mtor_geno_plot <- ggplot(ready1,aes(x=Coefficient,y=Estimate))
mtor_geno_plot <- mtor_geno_plot + geom_errorbar(aes(ymin=LowCI,ymax=HighCI),width=0.2) + geom_point(shape=16,size=3,fill="white") + geom_hline(yintercept=0,linetype="dashed") + scale_fill_discrete(name="") + xlab("") + ylab("Estimate") + coord_flip() + theme_classic() + theme(text=element_text(size=12)) + ggtitle("MTOR")

pdf("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Figures/genotyping_figure_strict/mtor_genotyping_eff_modeling.pdf", width = 10, height =4)
print(mtor_geno_plot)
dev.off()
