# differential gene expression between pulp and perio in combined dataset
# ex cellTypeChar = "Epithelial"
diff_exp_by_subpop <- function( cellTypeChar ){
library(ggplot2)
# library(cowplot)
# theme_set(theme_cowplot())
Idents(object = merged_harmony) <- "groups_bysize"
cellType <- subset(merged_harmony, idents = cellTypeChar)
Idents(cellType) <- "condition"
# avg.cellType <- log1p(AverageExpression(cellType, verbose = FALSE, return.seurat = TRUE)$RNA)
# avg.cellType$gene <- rownames(avg.cellType)

# DoHeatmap(avg.cellType, features = unlist(TopFeatures(pbmc[["pca"]], balanced = TRUE)), size = 3, 
#     draw.lines = FALSE)

# genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
# p1 <- ggplot(avg.cellType, aes(perio, pulp)) + geom_point() + ggtitle("Epithelial cells")
# p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)

# Combine cell type with condition label
merged_harmony$celltype.cond = paste(Idents(merged_harmony), merged_harmony$condition, sep = "_")
DefaultAssay(merged_harmony) = "RNA"
Idents(merged_harmony) 	 	 = "celltype.cond"
epithelial_by_cond 		     = FindMarkers(merged_harmony, ident.1 = paste0(cellTypeChar, "_perio"), ident.2 = paste0(cellTypeChar, "_pulp"), verbose = FALSE)
head(epithelial_by_cond, n = 15)
epithelial_by_cond$gene 	 = rownames(epithelial_by_cond)
top_epithelial				 = top_n(epithelial_by_cond[order(epithelial_by_cond$p_val_adj),], n = 50, wt = avg_logFC)
# top_epithelial$gene

# cluster.averages <- AverageExpression(object = cellType, return.seurat = TRUE)

g = DoHeatmap(object = cellType, angle = 45, slot = "scale.data", features = top_epithelial$gene, group.by= "condition", group.colors = c("blue", "red") )
ggsave(plot = g, filename = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses", 
	paste0("Heatmap_", cellTypeChar,"_bycondition.pdf")), width=16, height=20)



top_epithelial

}

# After running the above analysis for epithelial, immune, MSCs and fibroblasts, save the results of the differential gene expression for each of the those cases:
# save(list=c("epi_perio_pulp", "immu_perio_pulp", "msc_perio_pulp", "fibro_perio_pulp"), file=file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses", "dge.RData"))



