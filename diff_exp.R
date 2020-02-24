# differential gene expression between pulp and perio in combined dataset
# ex cellTypeChar = "Epithelial"
library(ggplot2)

Idents(object = merged_harmony) <- "groups_bysize"
merged_harmony$celltype.cond = paste(Idents(merged_harmony), merged_harmony$condition, sep = "_")

diff_exp_by_subpop <- function( cellTypeChar ){

cellType <- subset(merged_harmony, idents = cellTypeChar)

DefaultAssay(cellType) 		 = "RNA"
Idents(cellType) 	 	 	 = "celltype.cond"
cellType_by_cond 		     = FindMarkers(cellType, ident.1 = paste0(cellTypeChar, "_perio"), ident.2 = paste0(cellTypeChar, "_pulp"), verbose = FALSE)
head(cellType_by_cond, n = 15)
cellType_by_cond$gene 	 	 = rownames(cellType_by_cond)
top_cellType1				 = top_n(cellType_by_cond[order(cellType_by_cond$p_val_adj),], n = 50, wt = avg_logFC) # overexpressed in ident.1 (perio)
top_cellType2				 = top_n(cellType_by_cond[order(cellType_by_cond$p_val_adj),], n = -50, wt = avg_logFC) # overexp. in ident.2 (pulp) 

# test = rbind.data.frame(top_cellType1, top_cellType2)

g = DoHeatmap(object = cellType, angle = 45, slot = "scale.data", features = c(top_cellType1$gene,top_cellType2$gene), group.by= "condition", group.colors = c("blue", "red") )
	ggsave(plot = g, filename = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses", 
	paste0("Heatmap_", cellTypeChar,"_bycondition.pdf")), width=16, height=20)

list(perio_high = top_cellType1, pulp_high = top_cellType2)

}

# To run:
# fibro_perio_pulp 		= diff_exp_by_subpop("Fibroblasts")
# epi_perio_pulp 		= diff_exp_by_subpop("Epithelial")
# msc_perio_pulp 		= diff_exp_by_subpop("MSC")
# immu_perio_pulp 		= diff_exp_by_subpop("Immune")

# After running the above analysis for epithelial, immune, MSCs and fibroblasts, save the results of the differential gene expression for each of the those cases:
# save(list=c("epi_perio_pulp", "immu_perio_pulp", "msc_perio_pulp", "fibro_perio_pulp"), file=file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses", "dge.RData"))



