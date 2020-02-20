# UMAPs in supplement
library(ggplot2)
library(Seurat)

#####################################################################################################################################################
############################################################ Perio ##################################################################################
#####################################################################################################################################################

CTasgns_perio 				= read.csv("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPerio/PerioCTdefs.csv", header = FALSE)
names(CTasgns_perio)		= c("cluster", 'metasubtype','cellType', 'colors')
Idents(object = all_perio)  = "groups_bysize"

################################################### Dimplots and Feature plots by Cell Type ################################################################
cell_types = c("Fibroblasts", "Endothelial", "Immune", "Epithelial")
genes_lists = list("Fibroblasts" = c("COL1A1", "OMD", "OGN", "PDGFRA"),"Endothelial" = c("PECAM", "EMCM","ACKR1", "IGFBP3", "CLDN5"), "Immune" = c("PTPRC",
"CSF1R","NKG7","GZMA","GZMH","CCL5", "LTB","IL7R","CD3E","CD4","CD8"),  "Epithelial" = c("FDCSP","S100A8","KRT14","SLPI","TUBA1B","ACTB","IL1B","IL1A"))  # Gene lists provided by Pierfrancesco

for (ss in 1:length(cell_types)){

ll 							= cell_types[ss]
subtype 					= subset(all_perio, idents = ll)
subtype_ids 				= which(CTasgns_perio$metaCellType==ll)
subtype_labels 				= CTasgns_perio$cellType[subtype_ids]
subtype_colors				= CTasgns_perio$colors[subtype_ids]

	DimPlot(subtype, reduction = "umap", label = TRUE, pt.size = 3, group.by = "ordered_clustering") +
	scale_colour_manual('Clusters', values = as.character(subtype_colors), labels = subtype_labels)
	ggsave(filename = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPerio", paste0(ll ,".pdf")), width=16, height=20)

	# Marker expressions
	FeaturePlot(subtype,features = genes_lists[[ll]], 
		reduction = "umap", order = T)
	ggsave(filename = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPerio", paste0(ll,"_Features.pdf")), width=20, height=20)

} 

#####################################################################################################################################################
############################################################# Pulp ##################################################################################
#####################################################################################################################################################

CTasgns_pulp 				= read.csv("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPulp/PulpCTdefs.csv", header = FALSE)
names(CTasgns_pulp)		    = c("cluster", 'metasubtype','cellType', 'colors')
Idents(object = all_healthy)  = "groups_bysize"

################################################### Dimplots and Feature plots by Cell Type #############################################################
cell_types = c("Fibroblasts", "Endothelial", "Immune", "Epithelial", "Myelinating ScCs", "Non-Myelinating ScCs")
genes_lists_pulp = list("Fibroblasts" = c("COL1A1","IGFB5", "CXCL14","PTN","OMD","COCH","GOLIM4","CTNNB1"),"Endothelial" = c("PECAM", "EMCM","ACKR1", "RGCC","INSR","IGFBP3", "CLDN5", "EDN1", "FABP4", "POSTN"),
	"Immune" = c("PTPRC","CSF1R","NKG7","GZMA","GZMH","CCL5", "LTB","IL7R","CD3E","CD4","CD8"),  "Epithelial" = c("FDCSP","S100A8","S100A9", "KRT14","SFN"), 
	"Myelinating ScCs" = c("S100B","MBP","MPZ","PRX","GFRA3","PLP1","ANXA1"), "Non-Myelinating ScCs" = c("S100B","MBP","MPZ","PRX","GFRA3","PLP1","ANXA1"))  # Gene lists provided by Pierfrancesco

for (ss in 1:length(cell_types)){

ll 							= cell_types[ss]
subtype 					= subset(all_healthy, idents = ll)
subtype_ids 				= which(CTasgns_pulp$metaCellType==ll)
subtype_labels 				= CTasgns_pulp$cellType[subtype_ids]
subtype_colors				= CTasgns_pulp$colors[subtype_ids]

	DimPlot(subtype, reduction = "umap", label = TRUE, pt.size = 3, group.by = "ordered_clustering") +
	scale_colour_manual('Clusters', values = as.character(subtype_colors), labels = subtype_labels)
	ggsave(filename = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPulp", paste0(ll ,".pdf")), width=16, height=20)

	# Marker expressions
	FeaturePlot(subtype,features = genes_lists_pulp[[ll]], 
		reduction = "umap", order = T)
	ggsave(filename = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPulp", paste0(ll,"_Features.pdf")), width=10, height=10)

} 


