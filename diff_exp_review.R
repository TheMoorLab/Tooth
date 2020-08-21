# differential gene expression between pulp and perio in combined dataset
# ex cellTypeChar = "Epithelial"
library(ggplot2)
library(tidyverse)


# load merged_harmony seurat object
load("~/Downloads/merged2.rds")

# create cell type by condition meta.data
Idents(object = merged_harmony) <- "groups_bysize"
merged_harmony$celltype.cond 	= paste(Idents(merged_harmony), merged_harmony$condition, sep = "_")


# Find markers that are differentially expressed between Perio and Pulp MSCs
DefaultAssay(merged_harmony) 		= "RNA"
Idents(merged_harmony) 	 	 	 	= "celltype.cond"

msc_perio_pulp 			= FindMarkers(object = merged_harmony, ident.1 = 'MSC_perio', ident.2 = 'MSC_pulp', min.pct = 0.10) %>% rownames_to_column("gene")

geneset_msc_DE			= msc_perio_pulp %>% filter(p_val_adj <= 0.05 & abs(avg_logFC) >= 0.25) %>% pull(gene)

save.image(file = "merged202008.Rdata")

