########################
#I/O
########################
library(tidyverse)
library(ggplot2)
library(data.table)
library(ggExtra)

sample_names <- c("e174", "FC5501", "FC5801")

#read in genotyping data 
loose_data_list <- list()
loose_data_list[["e174"]] <- fread("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Genotyping/e174/scMosaicHunter_results_e174_using_collapsed_10X_umi_with_custom_VAF_looseGenotyping.txt", header = TRUE, sep = "\t")
loose_data_list[["FC5501"]] <- fread("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Genotyping/FC5501/scMosaicHunter_results_FC5501_using_collapsed_10X_umi_with_custom_VAF_looseGenotyping.txt", header = TRUE, sep = "\t")
loose_data_list[["FC5801"]] <- fread("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Genotyping/FC5801/scMosaicHunter_results_FC5801_using_collapsed_10X_umi_with_custom_VAF_looseGenotyping.txt", header = TRUE, sep = "\t")

loose_data_input_list <- list()
loose_data_input_list[["e174"]] <- fread("/n/groups/walsh/indData/Maya/FCD_project/raw_data/Edward_GOTEN_Data/e174.all.csv", header = TRUE, sep = ",") %>% as.data.frame()
loose_data_input_list[["FC5501"]] <- fread("/n/groups/walsh/indData/Maya/FCD_project/raw_data/Edward_GOTEN_Data/FC5501.all.csv", header = TRUE, sep = ",") %>% as.data.frame
loose_data_input_list[["FC5801"]] <- fread("/n/groups/walsh/indData/Maya/FCD_project/raw_data/Edward_GOTEN_Data/FC5801.all.csv", header = TRUE, sep = ",") %>% as.data.frame()
loose_data_input_list_filtered <- lapply(sample_names, function(x) loose_data_input_list[[x]] %>% filter(best.bc %in% (loose_data_list[[x]])$best.bc))

strict_data_list <- list()
strict_data_list[["e174"]] <- fread("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Genotyping/e174/scMosaicHunter_results_e174_using_collapsed_10X_umi_with_custom_VAF_strictGenotyping.txt", header = TRUE, sep = "\t")
strict_data_list[["FC5501"]] <- fread("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Genotyping/FC5501/scMosaicHunter_results_FC5501_using_collapsed_10X_umi_with_custom_VAF_strictGenotyping.txt", header = TRUE, sep = "\t")
strict_data_list[["FC5801"]] <- fread("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Genotyping/FC5801/scMosaicHunter_results_FC5801_using_collapsed_10X_umi_with_custom_VAF_strictGenotyping.txt", header = TRUE, sep = "\t")

strict_data_input_list <- list()
strict_data_input_list[["e174"]] <- fread("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Genotyping/e174/e174_regenerated_final_genotype_list_based_on_local_alignment.txt", header = TRUE) %>% as.data.frame()
strict_data_input_list[["FC5501"]] <- fread("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Genotyping/FC5501/FC5501_regenerated_final_genotype_list_based_on_local_alignment.txt", header = TRUE) %>% as.data.frame
strict_data_input_list[["FC5801"]] <- fread("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Genotyping/FC5801/FC5801_regenerated_final_genotype_list_based_on_local_alignment.txt", header = TRUE) %>% as.data.frame()
strict_data_input_list_filtered <- lapply(sample_names, function(x) strict_data_input_list[[x]] %>% filter(best.bc %in% (strict_data_list[[x]])$best.bc))

#read in genotyping metadata
loose_metadata <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Seurat_Objects/9_processed_seurat_obj_integrated_azimuthAnnotations_ourAnnotations_doubletsFiltered_FCD1Filtered_looseGenotyping_metadata.RDS")
strict_metadata <- readRDS("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Seurat_Objects/9_processed_seurat_obj_integrated_azimuthAnnotations_ourAnnotations_doubletsFiltered_FCD1Filtered_strictGenotyping_metadata.RDS")

########################
#LOOSE GENOTYPING DATA
########################
#add in number of umis per nucleus and mrf 
loose_data_list <- lapply(loose_data_list, function(x) x %>% 
                                                        mutate(numRefUMIs = nchar(ref.bases)) %>% 
                                                        mutate(numAltUMIs = nchar(alt.bases)) %>% 
                                                        mutate(MRF = numAltUMIs/(numAltUMIs + numRefUMIs)) %>% 
                                                        dplyr::select(best.bc, genotype, genoqual, numRefUMIs, numAltUMIs, MRF) %>% 
                                                        as.data.frame() %>% 
                                                        mutate(totalUMIs = numRefUMIs + numAltUMIs))

#plot
loose_data_list <- lapply(seq_along(loose_data_list), function(x) loose_data_list[[x]] %>% mutate(sample = names(loose_data_list)[x]))
loose_data_df <- do.call(rbind, loose_data_list)

base_plot <- ggplot(loose_data_df, aes(x = totalUMIs, y = MRF, color = genotype)) +
  geom_point(size = 3, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)) + 
      scale_color_manual(
      values = c("hetero" = "#FF0000", "ref-hom" = "#FBB03B")
    ) +
    labs(
      x = "Total Number of UMIs",
      y = "Mutant Read Fraction",
      color = "Genotype",
      title = "High Sensitivity Genotyping Approach"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    )

 final_plot <- ggMarginal(
    base_plot,
    type = "histogram",
    groupFill = TRUE,
    margins = "both"
  )
  
  pdf(paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Genotyping/loose_genotyping_MRF.pdf"))
  print(base_plot)
  dev.off()

#summary statistics 
loose_data_df %>% group_by(genotype) %>% summarize(med_MRF = median(MRF), med_numUMIs = median(totalUMIs))
loose_data_input_list_filtered <- lapply(seq_along(loose_data_input_list_filtered), function(x) loose_data_input_list_filtered[[x]] %>% dplyr::select(best.bc, estimated.umi) %>% mutate(sample = sample_names[x]))
loose_data_input_list_filtered <- lapply(seq_along(loose_data_input_list_filtered), function(x) merge(loose_data_input_list_filtered[[x]], loose_data_list[[x]], by = "best.bc"))
loose_data_input_df <- do.call(rbind, loose_data_input_list_filtered)
loose_data_input_df %>% group_by(best.bc, estimated.umi, genotype) %>% summarize(n = n()) %>% ungroup() %>% group_by(genotype) %>% summarize(medNumReadsPerUMI = median(n))
table((loose_data_input_df %>% filter(genotype == "ref-hom") %>% dplyr::select(best.bc, totalUMIs) %>% unique)$totalUMIs)
table((loose_data_input_df %>% filter(genotype == "ref-hom") %>% dplyr::select(best.bc, totalUMIs) %>% unique)$totalUMIs > 2)
table((loose_data_input_df %>% filter(genotype == "hetero") %>% dplyr::select(best.bc, totalUMIs) %>% unique)$totalUMIs)
table((loose_data_input_df %>% filter(genotype == "hetero") %>% dplyr::select(best.bc, totalUMIs) %>% unique)$totalUMIs > 2)

table((loose_data_input_df %>% filter(genotype == "ref-hom") %>% group_by(best.bc, estimated.umi, genotype) %>% summarize(n = n()) %>% ungroup() %>% mutate(bc_umi = paste0(best.bc, "_", estimated.umi)))$n)
table((loose_data_input_df %>% filter(genotype == "ref-hom") %>% group_by(best.bc, estimated.umi, genotype) %>% summarize(n = n()) %>% ungroup() %>% mutate(bc_umi = paste0(best.bc, "_", estimated.umi)))$n > 2)
table((loose_data_input_df %>% filter(genotype == "hetero") %>% group_by(best.bc, estimated.umi, genotype) %>% summarize(n = n()) %>% ungroup() %>% mutate(bc_umi = paste0(best.bc, "_", estimated.umi)))$n)
table((loose_data_input_df %>% filter(genotype == "hetero") %>% group_by(best.bc, estimated.umi, genotype) %>% summarize(n = n()) %>% ungroup() %>% mutate(bc_umi = paste0(best.bc, "_", estimated.umi)))$n > 2)


########################
#STRICT GENOTYPING DATA
########################
#add in number of umis per nucleus and mrf 
strict_data_list <- lapply(strict_data_list, function(x) x %>% 
                                                        mutate(numRefUMIs = nchar(ref.bases)) %>% 
                                                        mutate(numAltUMIs = nchar(alt.bases)) %>% 
                                                        mutate(MRF = numAltUMIs/(numAltUMIs + numRefUMIs)) %>% 
                                                        dplyr::select(best.bc, genotype, genoqual, numRefUMIs, numAltUMIs, MRF) %>% 
                                                        as.data.frame() %>% 
                                                        mutate(totalUMIs = numRefUMIs + numAltUMIs))

#plot
strict_data_list <- lapply(seq_along(strict_data_list), function(x) strict_data_list[[x]] %>% mutate(sample = names(strict_data_list)[x]))
strict_data_df <- do.call(rbind, strict_data_list)

base_plot <- ggplot(strict_data_df, aes(x = totalUMIs, y = MRF, color = genotype)) +
  geom_point(size = 3, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)) + 
      scale_color_manual(
      values = c("hetero" = "#FF0000", "ref-hom" = "#FBB03B")
    ) +
    labs(
      x = "Total Number of UMIs",
      y = "Mutant Read Fraction",
      color = "Genotype",
      title = "High Specificity Genotyping Approach"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    )

 final_plot <- ggMarginal(
    base_plot,
    type = "histogram",
    groupFill = TRUE,
    margins = "both"
  )
  
  pdf(paste0("/n/groups/walsh/indData/Maya/FCD_project/analysis/2_Analyze_Full_Object/Output/Genotyping/strict_genotyping_MRF.pdf"))
  print(base_plot)
  dev.off()

#summary statistics 
strict_data_df %>% group_by(genotype) %>% summarize(med_MRF = median(MRF), med_numUMIs = median(totalUMIs))
strict_data_input_list_filtered <- lapply(seq_along(strict_data_input_list_filtered), function(x) strict_data_input_list_filtered[[x]] %>% dplyr::select(best.bc, estimated.umi) %>% mutate(sample = sample_names[x]))
strict_data_input_list_filtered <- lapply(seq_along(strict_data_input_list_filtered), function(x) merge(strict_data_input_list_filtered[[x]], strict_data_list[[x]], by = "best.bc"))
strict_data_input_df <- do.call(rbind, strict_data_input_list_filtered)
strict_data_input_df %>% group_by(best.bc, estimated.umi, genotype) %>% summarize(n = n()) %>% ungroup() %>% group_by(genotype) %>% summarize(medNumReadsPerUMI = median(n))
table((strict_data_input_df %>% filter(genotype == "ref-hom") %>% dplyr::select(best.bc, totalUMIs) %>% unique)$totalUMIs)
table((strict_data_input_df %>% filter(genotype == "ref-hom") %>% dplyr::select(best.bc, totalUMIs) %>% unique)$totalUMIs > 2)
table((strict_data_input_df %>% filter(genotype == "hetero") %>% dplyr::select(best.bc, totalUMIs) %>% unique)$totalUMIs)
table((strict_data_input_df %>% filter(genotype == "hetero") %>% dplyr::select(best.bc, totalUMIs) %>% unique)$totalUMIs > 2)

table((strict_data_input_df %>% filter(genotype == "ref-hom") %>% group_by(best.bc, estimated.umi, genotype) %>% summarize(n = n()) %>% ungroup() %>% mutate(bc_umi = paste0(best.bc, "_", estimated.umi)))$n)
table((strict_data_input_df %>% filter(genotype == "ref-hom") %>% group_by(best.bc, estimated.umi, genotype) %>% summarize(n = n()) %>% ungroup() %>% mutate(bc_umi = paste0(best.bc, "_", estimated.umi)))$n > 2)
table((strict_data_input_df %>% filter(genotype == "hetero") %>% group_by(best.bc, estimated.umi, genotype) %>% summarize(n = n()) %>% ungroup() %>% mutate(bc_umi = paste0(best.bc, "_", estimated.umi)))$n)
table((strict_data_input_df %>% filter(genotype == "hetero") %>% group_by(best.bc, estimated.umi, genotype) %>% summarize(n = n()) %>% ungroup() %>% mutate(bc_umi = paste0(best.bc, "_", estimated.umi)))$n > 2)
