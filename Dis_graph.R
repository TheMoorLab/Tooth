# Jaccard distance plotted tih force-directed graph layout


CTasgns 			= read.csv("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/CellTypesDefinitions.csv", header = FALSE)
names(CTasgns)		= c("cluster", 'metaCellType','cellType', 'colors')

Idents(object = merged_harmony) <- "groups_bysize"
av.ex_top 			= AverageExpression(merged_harmony, assay= "RNA", add.ident = 'condition', features =  merged_harmony@assays$SCT@var.features)
# av.ex.SCT2 			= t(av.ex_top)

cluster_dist2   		= dist(t(av.ex_top) , method = "ejaccard")
dist_mat 		    = as.matrix(cluster_dist)
