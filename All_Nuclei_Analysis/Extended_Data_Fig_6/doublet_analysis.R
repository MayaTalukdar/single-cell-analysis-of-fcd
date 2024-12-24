library(ggrepel) 
library(ggplot2)
library(dplyr)
library(cowplot)

###################
#doublet analysis
###################
merged.Seurat.obj <- subset(merged.Seurat.obj, subset=origin=='OurData')
merged.Seurat.obj <- subset(merged.Seurat.obj,subset=status=='case')
merged.Seurat.obj@meta.data$is_doublet_class <- ifelse(
  merged.Seurat.obj@meta.data$is_putative_doublet == 'TRUE',
  "doublet",
  "singlet"
)

#examine putative doublets by cell type
cowplot::plot_grid(ncol = 2, DimPlot(merged.Seurat.obj, reduction = "umap.harmony", group.by = "is_putative_doublet"),
                              DimPlot(merged.Seurat.obj, reduction = "umap.harmony", group.by = "predicted.subclass"))

#examine putative doublets by samples                            
DimPlot(merged.Seurat.obj, reduction = "umap.harmony", group.by="is_doublet_class",split.by = "cleaned_ident")

#correlation analysis
dr <- as.data.frame(prop.table(table(merged.Seurat.obj$is_doublet_class,merged.Seurat.obj$cleaned_ident),2)["doublet",])
colnames(dr) <- "doublet_ratio"
cn <- as.data.frame(table(merged.Seurat.obj$cleaned_ident))
dr$cn <- cn$Freq
sinfo <- merged.Seurat.obj@meta.data[,c("cleaned_ident","origin")]%>%
  distinct()
sinfo <- sinfo[match(rownames(dr), sinfo$cleaned_ident), ]
dr$origin <- sinfo$origin
cor_result <- cor.test(dr$cn, dr$doublet_ratio, method = "pearson")
r_value <- cor_result$estimate
p_value <- cor_result$p.value

#correlation plot between doublet ratio and total cell number
ggplot(dr, aes(x = cn, y = doublet_ratio)) +
  geom_point(aes(color = origin)) + 
  geom_smooth(method = "lm", color = "grey", se = FALSE, linetype = "dashed") + 
  geom_text_repel(aes(label = rownames(dr)), size = 3) + 
  theme_bw() + 
  theme(
    aspect.ratio = 1,
    legend.title = element_blank(),
    legend.position = 'right'
  ) +
  coord_fixed() +
  ylab('doublet_ratio') +
  xlab('total cell number') +
  labs(
    title = paste0("coef = ", round(r_value, 2), 
                   ", p = ", signif(p_value, 3))
  )

vaf <- read.csv("vaf.csv", header = TRUE)
dr$id <- rownames(dr)
vaf <- merge(vaf,dr,by='id')
cor_result <- cor.test(vaf$vaf, vaf$doublet_ratio, method = "pearson")
r_value <- cor_result$estimate
p_value <- cor_result$p.value

#correlation plot between doublet ratio and vaf
ggplot(vaf, aes(x = vaf, y = doublet_ratio)) +
  geom_point(aes(color = origin)) +  
  geom_smooth(method = "lm", color = "grey", se = FALSE, linetype = "dashed") + 
  geom_text_repel(aes(label = id), size = 3) +  
  theme_bw() + 
  theme(
    aspect.ratio = 1,
    legend.title = element_blank(),
    legend.position = 'right'
  ) +
  coord_fixed() +
  ylab('doublet_ratio') +
  xlab('VAF') +
  labs(
    title = paste0("coef = ", round(r_value, 2), 
                   ", p = ", signif(p_value, 3))
  )


