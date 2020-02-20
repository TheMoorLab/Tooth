# For all perio and all pulp
# Umap of clusters
# heatmaps 
library(dplyr)
library(ggplot2)

###########################################################################################################################
################################################### Combined ##############################################################
###########################################################################################################################

load("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/r_output/merged_harmony_markers_20191204.rds")

Idents(object = merged_harmony) <- "clustering"

# Import clusters' cell type assignments and pre-defined custom color squeme
CTasgns_merged 				= read.csv("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/CellTypesDefinitions.csv", header = FALSE)
names(CTasgns_merged)		= c("cluster", 'metaCellType','cellType', 'colors')
new.cluster.ids 			= as.character(CTasgns_merged$cellType)
custom_color_squeme  		= as.character(CTasgns_merged$colors)

# Create cell identity based on meta cell type groups 
# This is important to be able to order levels below according to most populous meta cell type 
match_cl_ct 	= cbind.data.frame(cluster = CTasgns_merged$cluster, metaCT = CTasgns_merged$metaCellType)
cell_cl_id 		= merged_harmony@meta.data$clustering 
metaCT 			= match_cl_ct[match(cell_cl_id, match_cl_ct[,"cluster"]),2]
merged_harmony$groups = metaCT
metaCT_size 	= factor(metaCT, levels = names(sort(table(metaCT), decreasing =TRUE)))
merged_harmony$groups_bysize = metaCT_size

### Re-level cell identities according to DECREASING POP SIZE order
ordered_lev_mct 			= factor(CTasgns_merged$metaCellType, levels = names(sort(table(metaCT), decreasing =TRUE)))
mct_ids 					= order(ordered_lev_mct)
my_mct_levels 				= levels(merged_harmony)[mct_ids]
merged_harmony@active.ident = factor(x = merged_harmony@active.ident, levels = my_mct_levels)

pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses", "Perio_Pulp_Clusters.pdf"), width=16, height=12)
par(xpd=TRUE)
	DimPlot(merged_harmony, reduction = "umap", label = TRUE, pt.size = 0.5)+
	scale_colour_manual('Clusters', values = custom_color_squeme[mct_ids], labels = new.cluster.ids[mct_ids])
dev.off()

# Clusters by perio/pulp identity
pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses", "Perio_Pulp_Clusters_bysampleID.pdf"), width=16, height=12)
par(xpd=TRUE)
	DimPlot(merged_harmony, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = "orig.ident")
	# +
	# scale_colour_manual('Clusters', values = custom_color_squeme[mct_ids], labels = new.cluster.ids[mct_ids])
dev.off()

# By perio/pulp only (combining all patient from each condition)
pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses", "Perio_Pulp_Clusters_bycondition.pdf"), width=16, height=12)
par(xpd=TRUE)
	DimPlot(merged_harmony, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = "condition")
dev.off()


##########################################################################################################################
################################################### Perio ################################################################
##########################################################################################################################

load('/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/r_output/all_perio_markers_20191204.rds')


pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPerio", "Clusters.pdf"))
par(xpd=TRUE)
		DimPlot(all_perio, reduction = "umap",
        group.by = "clustering", pt.size = 0.5, label = FALSE, repel = TRUE) 
dev.off()


pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPerio", "Features.pdf"), width=16, height=10)
par(xpd=TRUE)
	FeaturePlot(all_perio,features=c("THY1","MYH11","NOTCH3","COL1A1","KRT14","PECAM1","PTPRC","MBP","GFRA3","MS4A1","SPOCK3", "HBB", "ADGRE1","SOX10"), 
		reduction = "umap",order = T)
dev.off()

# Clusters named according to cell type
# IMPORTANT: Ensure cell identities are clusters before running the code below
Idents(object = all_perio) <- "clustering"

# Import clusters' cell type assignments and pre-defined custom color squeme
CTasgns_perio 				= read.csv("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPerio/PerioCTdefs.csv", header = FALSE)
names(CTasgns_perio)		= c("cluster", 'metaCellType','cellType', 'colors')
new.cluster.ids 			= as.character(CTasgns_perio$cellType)
custom_color_squeme  		= as.character(CTasgns_perio$colors)

# create cell identity based on meta cell type groups with match. 
# This is important to be able to order levels below according to most populous meta cell type 
match_cl_ct 	= cbind.data.frame(cluster = CTasgns_perio$cluster, metaCT = CTasgns_perio$metaCellType)
cell_cl_id 		= all_perio@meta.data$clustering 
metaCT 			= match_cl_ct[match(cell_cl_id, match_cl_ct[,"cluster"]),2]
all_perio$groups = metaCT
metaCT_size 	= factor(metaCT, levels = names(sort(table(metaCT), decreasing =TRUE)))
all_perio$groups_bysize = metaCT_size

### Re-level cell identities according to DECREASING POP SIZE order
ordered_lev_mct 			= factor(CTasgns_perio$metaCellType, levels = names(sort(table(metaCT), decreasing =TRUE)))
mct_ids 					= order(ordered_lev_mct)
my_mct_levels 				= levels(all_perio)[mct_ids]
all_perio@active.ident 		= factor(x = all_perio@active.ident, levels = my_mct_levels)

#########################################################################################################################################################
# PLOT CLUSTERS WITH LABELS ACCORDING TO DECREASING POP SIZE order
# Cluster labels
pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPerio", "Clusters_ct_labels_pop_size_ordered_pt1_5.pdf"), width=12, height=12)
par(xpd=TRUE)
	DimPlot(all_perio, reduction = "umap", label = TRUE, pt.size = 1.5)+
	scale_colour_manual('Clusters', values = custom_color_squeme[mct_ids], labels = new.cluster.ids[mct_ids])
dev.off()

# Cluster numbers 
pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPerio", "Clusters_ct_nums_pop_size_ordered_new.pdf"), width=12, height=12)
par(xpd=TRUE)
	DimPlot(all_perio, reduction = "umap", label = TRUE, pt.size = 1)+
	scale_colour_manual('Clusters', values = custom_color_squeme[mct_ids])
dev.off()

# Create an array with the subclusters to add to meta.data (useful for future plots)
match_cl_ct 	= cbind.data.frame(cluster = CTasgns_perio$cluster, metaCT = CTasgns_perio$cellType)
cell_cl_id 		= all_perio@meta.data$clustering 
metaCT 			= match_cl_ct[match(cell_cl_id, match_cl_ct[,"cluster"]),2]


#########################################################################################################################################################
# HEATMAPS
# For correct heatmap, FindAllMarkers has to be applied to the Seurat object with the correct RE-LEVELED cluster groups
all.markers.perio = FindAllMarkers(all_perio, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10perio 		  = all.markers.perio %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top5perio 		  = all.markers.perio %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPerio", "Heatmap_pop_size_ordered_genes.pdf"), width=24, height=12)
par(xpd=TRUE)
DoHeatmap(object = all_perio, angle = 45, slot = "scale.data", features = top10perio$gene, group.by= "groups_bysize")
dev.off()

# Multi-bar heatmap: 
meta_color_squeme 						= c(custom_color_squeme[mct_ids][1], custom_color_squeme[mct_ids][6], custom_color_squeme[mct_ids][8], custom_color_squeme[mct_ids][9], custom_color_squeme[mct_ids][10], custom_color_squeme[mct_ids][13], custom_color_squeme[mct_ids][14])
cols.use 								= list(ordered_clustering=custom_color_squeme[mct_ids], groups_bysize = meta_color_squeme )
all_perio@meta.data$ordered_clustering  = all_perio@active.ident

pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPerio", "Heatmap_pop_size_ordered_genes_double_bar.pdf"), width=24, height=28)
par(xpd=TRUE)
DoMultiBarHeatmap(all_perio, features=top5perio$gene , group.by='groups_bysize', additional.group.by = 'ordered_clustering', additional.group.sort.by = 'ordered_clustering', cols.use=cols.use ) +
  theme(text = element_text(size = 20))+ guides( color = FALSE, size = FALSE)
  	# , legend.position = "none") 
dev.off()


pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPerio", "Heatmap_pop_size_ordered_subcellclusters.pdf"), width=24, height=12)
par(xpd=TRUE)
DoHeatmap(object = all_perio, angle = 45, slot = "scale.data", features = top10perio$gene, group.by= "ordered_clustering")
dev.off()

##########################################################################################################################
#################################################### Pulp ################################################################
##########################################################################################################################

load('/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/r_output/all_healthy_markers_20191129.rds')

# Clusters numbered
pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPulp", "Clusters.pdf"))
par(xpd=TRUE)
		DimPlot(all_healthy, reduction = "umap",
        group.by = "clustering", pt.size = 0.5, label = TRUE, repel = TRUE) 
dev.off()

# Marker expressions
pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPulp", "Features.pdf"), width=16, height=10)
par(xpd=TRUE)
	FeaturePlot(all_healthy,features=c("THY1","MYH11","COL1A1","KRT14","PECAM1","PTPRC","MBP","GFRA3","MS4A1","SPOCK3", "SOX10"), 
		reduction = "umap", order = T)
dev.off()

# Clusters named according to cell type
# Ensure cell identities are clusters before running the code below
Idents(object = all_healthy) <- "clustering"

# Import clusters' cell type assignments and pre-defined custom color squeme
CTasgns_pulp 				= read.csv("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPulp/PulpCTdefs.csv", header = FALSE)
names(CTasgns_pulp)			= c("cluster", 'metaCellType','cellType', 'colors')
new.cluster.ids 			= as.character(CTasgns_pulp$cellType)
custom_color_squeme  		= as.character(CTasgns_pulp$colors)

# Create cell identity based on meta cell type groups with match. 
# This is important to be able to order levels below according to most populous meta cell type 
match_cl_ct 	= cbind.data.frame(cluster = CTasgns_pulp$cluster, metaCT = CTasgns_pulp$metaCellType)
cell_cl_id 		= all_healthy@meta.data$clustering 
metaCT 			= match_cl_ct[match(cell_cl_id, match_cl_ct[,"cluster"]),2]
all_healthy$groups = metaCT
metaCT_size 	= factor(metaCT, levels = names(sort(table(metaCT), decreasing =TRUE)))
all_healthy$groups_bysize = metaCT_size

### Re-level cell identities according to DECREASING POP SIZE order
ordered_lev_mct 			= factor(CTasgns_pulp$metaCellType, levels = names(sort(table(metaCT), decreasing =TRUE)))
mct_ids 					= order(ordered_lev_mct)
my_mct_levels 				= levels(all_healthy)[mct_ids]
all_healthy@active.ident 	= factor(x = all_healthy@active.ident, levels = my_mct_levels)

# PLOT CLUSTERS WITH LABELS ACCORDING TO DECREASING POP SIZE order
# Cluster labels
pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPulp", "Clusters_ct_labels_pop_size_ordered_new.pdf"), width=12, height=12)
par(xpd=TRUE)
	DimPlot(all_healthy, reduction = "umap", label = TRUE, pt.size = 0.5)+
	scale_colour_manual('Clusters', values = custom_color_squeme[mct_ids], labels = new.cluster.ids[mct_ids])
dev.off()

# Cluster numbers 
pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPulp", "Clusters_ct_nums_pop_size_ordered_new.pdf"), width=12, height=12)
par(xpd=TRUE)
	DimPlot(all_healthy, reduction = "umap", label = TRUE, pt.size = 0.5)+
	scale_colour_manual('Clusters', values = custom_color_squeme[mct_ids],)
dev.off()

#########################################################################################################################################################
# HEATMAPS Ordered cells and genes
# For correct heatmap, FindAllMarkers has to be applied to the Seurat object with the correct RE-LEVELED cluster groups
all.markers.releveled = FindAllMarkers(all_healthy, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# top20 				  = all.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
DoHeatmap(object = all_healthy, genes.use = top20$gene, slim.col.label = TRUE, remove.key = TRUE)
top10 				  = all.markers.releveled %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

top5 				  = all.markers.releveled %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)


pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPulp", "Heatmap_pop_size_ordered_genes.pdf"), width=24, height=12)
par(xpd=TRUE)
DoHeatmap(object = all_healthy, angle = 45, slot = "scale.data", features = top10$gene, group.by= "groups_bysize")
dev.off()

meta_color_squeme_pulp 					  = c(custom_color_squeme[mct_ids][1], custom_color_squeme[mct_ids][7], custom_color_squeme[mct_ids][12], custom_color_squeme[mct_ids][15], custom_color_squeme[mct_ids][20], custom_color_squeme[mct_ids][21], custom_color_squeme[mct_ids][22], custom_color_squeme[mct_ids][23], custom_color_squeme[mct_ids][24])
cols.use 								  = list(ordered_clustering=custom_color_squeme[mct_ids], groups_bysize = meta_color_squeme_pulp )
all_healthy@meta.data$ordered_clustering  = all_healthy@active.ident

pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPulp", "Heatmap_pop_size_ordered_genes_double_bar.pdf"), width=24, height=28)
par(xpd=TRUE)
DoMultiBarHeatmap(all_healthy, features=top5$gene , size = 6, group.by='groups_bysize', additional.group.by = 'ordered_clustering', additional.group.sort.by = 'ordered_clustering', cols.use=cols.use ) +
  theme(text = element_text(size = 20))+ guides( color = FALSE, size = FALSE)
  	# , legend.position = "none") 
dev.off()

pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPulp", "Heatmap_pop_size_ordered_clusters.pdf"), width=24, height=20)
par(xpd=TRUE)
DoHeatmap(object = all_healthy, angle = 45, slot = "scale.data", features = top10$gene, group.by= "ordered_clustering")
dev.off()

################################################################## Save data #############################################################################

save.image(file = "perio_pulp_merged_20191220.Rdata")

##########################################################################################################################################################
############################################################## COMBINED PULP AND PERIO ###################################################################
# PROPORTIONS
# FOR MERGED DATASET OF PULP AND PERIO
# Create perio/pulp only annotation without distinguishing by sample Perio 1, 2, etc...
red_orig.ident 		 = gsub('[[:digit:]]+', '', merged_harmony@meta.data$orig.ident)
red_orig.ident2 	 = gsub('v', '', red_orig.odent)
meta.data$condition = red_orig.ident2

table(paste0(merged_harmony@meta.data$condition, merged_harmony@meta.data$groups_bysize))

#basic plot of clusters by replicate
ggplot(merged_harmony@meta.data, aes(x=groups_bysize, fill=orig.ident)) + 
	geom_bar(position = "fill")+
 	labs(x="Cell type", y="Proportion of sample",
       title="Proportion of sample for each cell type", fill = "Sample") + 
  theme(axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 1), axis.title=element_text(size=14,face="bold"))
ggsave(filename = "/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/Combined_proportions_samples.pdf")


#plot as proportion or percentage of cluster
ggplot(merged_harmony@meta.data, aes(x=groups_bysize, fill=condition)) + 
	geom_bar(position = "fill")+
 	labs(x="Cell type", y="Proportion of sample",
       title="Proportion of sample for each cell type", fill = "Sample") + 
  theme(axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 1), axis.title=element_text(size=14,face="bold"))
ggsave(filename = "/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/Combined_proportions_perio_pulp.pdf")

#plot as proportion or percentage of cluster normalized by number of cells in perio or pulp ( since it's disproportionate towards more pulp)
melted_celltype_condition  		   = cbind.data.frame(cellType =merged_harmony@meta.data$groups_bysize, condition=merged_harmony@meta.data$condition)
melted_celltype_by_condition 	   = as.data.frame.matrix(table(melted_celltype_condition))
cellType_normalized_by_condition   = as.data.frame(t(melted_celltype_by_condition)/rowSums(t(melted_celltype_by_condition)))

datm2 <- cellType_normalized_by_condition %>% 
  mutate(ind = factor(row_number())) %>% 
  gather(variable, value, -ind)

ggplot(datm2, aes(x = variable, y = value, fill = ind)) + 
    geom_bar(position = "fill",stat = "identity") +
    # or:
    # geom_bar(position = position_fill(), stat = "identity") 
    scale_y_continuous(labels = scales::percent_format())+
     	labs(x="Cell type", y="Normalized proportion",
       title="Proportion of cell type from sample (normalized)", fill = "Sample") +
     	  scale_fill_manual(labels = c("Perio", "Pulp"), values = c("blue", "red")) +
     	   theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/Normalized_proportions_perio_pulp.pdf")



##################################################################################################################################
# # sort cell types into corresponding groups (all endo together, all fibro, etc...)
# Idents(object = all_healthy) <- "clustering"
# sorted_ct 					= sort(new.cluster.ids, index.return=TRUE)
# my_levels 					= levels(all_healthy)[sorted_ct$ix]

# ### Re-level cell identities according to our pre-defined ALPHABETICAL order groups
# all_healthy@active.ident 	= factor(x = all_healthy@active.ident, levels = my_levels)

# # with cluster labels
# pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPulp", "Clusters_ct_labels.pdf"), width=12, height=12)
# par(xpd=TRUE)
# 	DimPlot(all_healthy, reduction = "umap", label = TRUE, pt.size = 0.5)+
# 	# , group.by = 'ident') # default, but can change to color points by specific features (clustering, cell cycle stage...)
# 	scale_colour_manual('Clusters', values = custom_color_squeme[sorted_ct$ix], labels = new.cluster.ids[sorted_ct$ix])
# 	# + NoLegend()
# dev.off()

# # With cluster numbers 
# pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPulp", "Clusters_ct_nums.pdf"), width=12, height=12)
# par(xpd=TRUE)
# 	DimPlot(all_healthy, reduction = "umap", label = TRUE, pt.size = 0.5)+
# 	# , group.by = 'ident') # default, but can change to color points by specific features (clustering, cell cycle stage...)
# 	scale_colour_manual('Clusters', values = custom_color_squeme[sorted_ct$ix])
# 	# + NoLegend()
# dev.off()
##################################################################################################################################



