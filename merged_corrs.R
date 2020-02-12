# Load merged data from Andreas:
library(corrplot)
library(amap) # distances
library(proxy) # more distances
library(qgraph) # force-directed graph layout
library(reshape2)


load("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/r_output/merged_harmony_markers_20191204.rds")

# Define cluster annotations:
pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses", "Perio_Pulp_celltypemarkers.pdf"),  width = 10, height = 8)
FeaturePlot(merged_harmony,features=c("THY1","MYH11","NOTCH3","COL1A1","KRT14","PECAM1","PTPRC","MBP","GFRA3","SPOCK3","MS4A1", "SOX10", "ADGRE1", "CD3E", "HBB"), reduction = "umap",order = T)
dev.off()

# Plot UMAP embedding of harmonized Perio and Pulp data clusters 
pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses", "Perio_Pulp_clustering.pdf"))
DimPlot(merged_harmony, label = TRUE, reduction = "umap", group.by = "clustering")
dev.off()

# # Pseudobulk corr
# av.ex     = AverageExpression(merged_harmony, add.ident = "orig.ident")
# av.ex.SCT = av.ex$SCT
# cor.exp   = as.data.frame(cor(av.ex.SCT))
# samples   = unique(gsub(".*_", "", row.names(cor.exp)))

# for (ii in 1:length(samples)){

# 	curr_samp 		 = samples[ii]
# 	curr_ids  		 = grep(curr_samp, row.names(cor.exp))
# 	curr_M 	  		 = cor.exp[curr_ids, curr_ids]

# 	clusters 	     = gsub("_.*","",row.names(curr_M))

# 	# Sort correlation matrix and its labels
# 	cids_ordered     = sort(as.numeric(clusters), index.return = TRUE)

# 	curr_M_sorted 	 = curr_M[cids_ordered$ix, cids_ordered$ix]
# 	row.names(curr_M_sorted)= paste("Cluster", cids_ordered$x)
# 	colnames(curr_M_sorted) = paste("Cluster", cids_ordered$x)

# 	pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/PseudobulkCorrs", paste0(curr_samp, ".pdf")))
# 	par(xpd=TRUE)
# 	corrplot(as.matrix(curr_M_sorted), method = "circle", type = "upper", 
# 		# order = "hclust", 
# 		title = curr_samp,
# 		mar = c(2, 2, 4, 4), number.cex = 0.5, 
# 		number.digits = 2, tl.col = 'black',
# 		tl.cex = 1)
# 	dev.off()
# }

# Pseudobulk corr on top 3000 variable features only
# av.ex_top = AverageExpression(merged_harmony, add.ident = "orig.ident", features =  merged_harmony@assays$SCT@var.features)
# av.ex.SCT = av.ex_top$SCT
# cor.exp   = as.data.frame(cor(av.ex.SCT))
# samples   = unique(gsub(".*_", "", row.names(cor.exp)))

# for (ii in 1:length(samples)){

# 	curr_samp 		 = samples[ii]
# 	curr_ids  		 = grep(curr_samp, row.names(cor.exp))
# 	curr_M 	  		 = cor.exp[curr_ids, curr_ids]

# 	clusters 	     = gsub("_.*","",row.names(curr_M))

# 	# Sort correlation matrix and its labels
# 	cids_ordered     = sort(as.numeric(clusters), index.return = TRUE)

# 	curr_M_sorted 	 = curr_M[cids_ordered$ix, cids_ordered$ix]
# 	row.names(curr_M_sorted)= paste("Cluster", cids_ordered$x)
# 	colnames(curr_M_sorted) = paste("Cluster", cids_ordered$x)

# 	pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/PseudobulkCorrs/TopFeats", paste0(curr_samp, ".pdf")))
# 	par(xpd=TRUE)
# 	corrplot(as.matrix(curr_M_sorted), method = "circle", type = "upper", 
# 		# order = "hclust", 
# 		title = curr_samp,
# 		mar = c(2, 2, 4, 4), number.cex = 0.5, 
# 		number.digits = 2, tl.col = 'black',
# 		tl.cex = 1)
# 	dev.off()
# }


# # Pseudobulk corr only on top variable features for AVERAGED all Perio and all Pulp 
# # Average expression on raw data
# CTasgns 			= read.csv("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/CellTypesDefinitions.csv", header = FALSE)
# names(CTasgns)		= c("cluster", 'cellType')

# # Normalize, scale and transform data
# merged_harmony  	= NormalizeData(merged_harmony, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
# merged_harmony 		= FindVariableFeatures(merged_harmony, assay = "RNA", selection.method = "vst", nfeatures = 1000)
# all.genes 			= rownames(merged_harmony)
# merged_harmony 		= ScaleData(merged_harmony, assay = "RNA", features = all.genes)

# # Make sure active identity is clustering" 
# Idents(merged_harmony)= "clustering"

# av.ex_top 			= AverageExpression(merged_harmony, assay = "RNA", add.ident = "orig.ident", features =  merged_harmony@assays$RNA@var.features)

# av.ex.SCT			= t(av.ex_top$RNA)

# conditions 			= c("perio", "pulp")
# conditions_array 	= character()

# cc_ids     			= lapply(conditions, function(a) grep(pattern = a, x = row.names(av.ex.SCT))) 

# # assign conditions to their respective locations in array
# conditions_array[cc_ids[[1]]] = conditions[1] 
# conditions_array[cc_ids[[2]]] = conditions[2]

# # Find cluster assignments
# clusters 			= gsub("_.*","",row.names(av.ex.SCT))

# # sub for imported cell type assignments
# CT_cluster  		= CTasgns[match(clusters, CTasgns[,"cluster"]),2]

# # Bind assignments and clusters to data 
# av_cl 				= cbind.data.frame(av.ex.SCT, cluster = CT_cluster, condition = conditions_array )

# # Aggregate data according to condition and cluster
# agg_cl 				= aggregate(.~cluster+condition, av_cl, mean)

# # Remove cluster and condition columns before running correlation
# data_agg 			=  agg_cl[-c(1,2)]

# # Estimate correlations across all subpopulations (all clusters of all conditions)
# cor.exp   			= cor(as.matrix(t(data_agg)))

# # Adjust correlation matrix labels
CT_cluster_condition= paste(agg_cl$condition, agg_cl$cluster)
# row.names(cor.exp) 	= CT_cluster_condition
# colnames(cor.exp) 	= CT_cluster_condition

# # Save correlation plots
# pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/PseudobulkCorrs", 
# 	"PerioPulpCorrelations_CT_top1k.pdf"))
# par(xpd=TRUE)
# corrplot(as.matrix(cor.exp), method = "circle", type = "upper", 
# 	order = "hclust", 
# 	mar = c(2, 2, 4, 4), number.cex = 0.5, 
# 	number.digits = 2, tl.col = 'black',
# 	tl.cex = 1)
# dev.off()


# Pseudobulk corr
# Top variable features for all Perio and all Pulp AVERAGED across patients in each group
# Average expression on SCT transformed data
# Import cluster assignments
CTasgns 			= read.csv("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/CellTypesDefinitions.csv", header = FALSE)
names(CTasgns)		= c("cluster", 'metaCellType','cellType', 'colors')

# av.ex_top 			= AverageExpression(merged_harmony, add.ident = "orig.ident", features =  merged_harmony@assays$SCT@var.features)
av.ex.SCT 			= t(av.ex_top$SCT)

conditions 			= c("perio", "pulp")
conditions_array 	= character()

cc_ids     			= lapply(conditions, function(a) grep(pattern = a, x = row.names(av.ex.SCT))) 

# assign conditions to their respective locations in array
conditions_array[cc_ids[[1]]] = conditions[1] 
conditions_array[cc_ids[[2]]] = conditions[2]

# Find cluster assignments
clusters 			= gsub("_.*","",row.names(av.ex.SCT))

# sub for imported cell type assignments
CT_cluster  		= CTasgns[match(clusters, CTasgns[,"cluster"]),2]

# Bind assignments and clusters to data 
av_cl 				= cbind.data.frame(av.ex.SCT, cluster = CT_cluster, condition = conditions_array )

# Aggregate data according to condition and cluster
agg_cl 				= aggregate(.~cluster+condition, av_cl, mean)

# Remove cluster and condition columns before running correlation
data_agg 			=  agg_cl[-c(1,2)]

# Aggregate cluster by cluster number
# av_cl_num 			= cbind.data.frame(av.ex.SCT, cluster = clusters, condition = conditions_array )
# agg_cl_num 			= aggregate(.~cluster+condition, av_cl_num, mean)

# # Estimate correlations across all subpopulations (all clusters of all conditions)
# cor.exp   			= cor(as.matrix(t(data_agg)))

# Estimate Distances between all subpopulations (all clusters of all conditions)
cluster_dist   		= dist(data_agg, method = "ejaccard")
dist_mat 		    = as.matrix(cluster_dist)

# Adjust correlation matrix labels 
# CT_cluster_condition= paste(agg_cl$condition, agg_cl$cluster)
# row.names(cor.exp) 	= CT_cluster_condition
# colnames(cor.exp) 	= CT_cluster_condition
CT_cluster_condition= paste(agg_cl$condition, agg_cl$cluster)


cc_collapsed_ids    = lapply(conditions, function(a) grep(pattern = a, x = CT_cluster_condition)) 
names(cc_collapsed_ids) = conditions

nodes_numbers 		= 1:length(CT_cluster_condition)
row.names(dist_mat) = nodes_numbers
colnames(dist_mat) 	= nodes_numbers

# Force directed graph layout
dist_mi <- 1-dist_mat  # 1 as qgraph takes similarity matrices as input
pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/PseudobulkCorrs", 
	"PerioPulpEjaccard_FDGL_mct_clusters2.pdf"), width=16, height=16)
# qgraph(dist_mi, layout='spring', vsize=3, groups = test, legend = TRUE)
qgraph(dist_mi, layout='spring', vsize=3, nodeNames = CT_cluster_condition, groups = cc_collapsed_ids, legend.mode="style2")
dev.off()

# HEATMAP visualization of the distance matrix dist_mat
row.names(dist_mi) 	= CT_cluster_condition
colnames(dist_mi) 	= CT_cluster_condition
melted_distmat 		= melt(dist_mi)

# ggplot(data = melted_distmat, aes(x=Var1, y=Var2, fill=value)) + geom_tile()

 g = ggplot(data = melted_distmat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0.5, limit = c(0,1), space = "Lab", 
   name="Similarity") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1), axis.title =element_blank())+
 coord_fixed()
ggsave(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/PseudobulkCorrs", 
	"Distance_mct_clusters2.pdf"), plot = g, width = 30, height = 20, units = "cm")

# Save correlation plots
# pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/PseudobulkCorrs", 
# 	"PerioPulpCorrelations_CT.pdf"))
# pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/PseudobulkCorrs", 
# 	"PerioPulpEjaccard_CT.pdf"))
# par(xpd=TRUE)
# # corrplot(as.matrix(cor.exp), method = "circle", type = "upper", 
# corrplot(as.matrix(dist_mat), method = "circle", type = "upper", 
# 	order = "hclust", 
# 	mar = c(2, 2, 4, 4), number.cex = 0.5, 
# 	number.digits = 2, tl.col = 'black',
# 	tl.cex = 1)
# dev.off()


# #######################
# Check pairwise relationships between clusters


# #######################
# Estimate correlations for each subset of the data matrix according to condition (perio or pulp)
# cc_sum_ids  = lapply(conditions, function(a) grep(pattern = a, x = agg_cl$condition)) 

# cor.exp   	= lapply(cc_sum_ids, function(a) cor(as.matrix(t(data_agg[a,]))))

# # assign correct cluster names to rownames 
# row.names(cor.exp[[1]]) 	= agg_cl$cluster[cc_sum_ids[[1]]]
# colnames(cor.exp[[1]]) 		= agg_cl$cluster[cc_sum_ids[[1]]]

# row.names(cor.exp[[2]]) 	= agg_cl$cluster[cc_sum_ids[[2]]]
# colnames(cor.exp[[2]]) 		= agg_cl$cluster[cc_sum_ids[[2]]]

# # Save correlation plots
# pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/PseudobulkCorrs", paste0(conditions[1],"_CT.pdf")))
# par(xpd=TRUE)
# corrplot(as.matrix(cor.exp[[1]]), method = "circle", type = "upper", 
# 	order = "hclust", 
# 	title = conditions[1],
# 	mar = c(2, 2, 4, 4), number.cex = 0.5, 
# 	number.digits = 2, tl.col = 'black',
# 	tl.cex = 1)
# dev.off()

# pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/PseudobulkCorrs", paste0(conditions[2],"_CT.pdf")))
# par(xpd=TRUE)
# corrplot(as.matrix(cor.exp[[2]]), method = "circle", type = "upper", 
# 	order = "hclust", 
# 	title = conditions[2],
# 	mar = c(2, 2, 4, 4), number.cex = 0.5, 
# 	number.digits = 2, tl.col = 'black',
# 	tl.cex = 1)
# dev.off()

# #######################

# # Bind assignments and clusters to data 
# av_cl 		= cbind.data.frame(av.ex.SCT, cluster = clusters, condition = conditions_array )

# # Aggregate data according to condition and cluster
# agg_cl 		= aggregate(.~cluster+condition, av_cl, mean)

# # Estimate correlations 
# cc_sum_ids  = lapply(conditions, function(a) grep(pattern = a, x = agg_cl$condition)) 
# data_agg 	=  agg_cl[-c(1,2)]
# cor.exp   	= lapply(cc_sum_ids, function(a) cor(as.matrix(t(data_agg[a,]))))

# # assign correct cluster names to rownames 
# row.names(cor.exp[[1]]) 	= agg_cl$cluster[cc_sum_ids[[1]]]
# colnames(cor.exp[[1]]) 		= agg_cl$cluster[cc_sum_ids[[1]]]

# row.names(cor.exp[[2]]) 	= agg_cl$cluster[cc_sum_ids[[2]]]
# colnames(cor.exp[[2]]) 		= agg_cl$cluster[cc_sum_ids[[2]]]

# # Save correlation plots
# pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/PseudobulkCorrs", "Pulp.pdf"))
# par(xpd=TRUE)
# corrplot(as.matrix(cor.exp[[1]]), method = "circle", type = "upper", 
# 	order = "hclust", 
# 	title = conditions[1],
# 	mar = c(2, 2, 4, 4), number.cex = 0.5, 
# 	number.digits = 2, tl.col = 'black',
# 	tl.cex = 1)
# dev.off()

# pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/PseudobulkCorrs", "Perio.pdf"))
# par(xpd=TRUE)
# corrplot(as.matrix(cor.exp[[2]]), method = "circle", type = "upper", 
# 	order = "hclust", 
# 	title = conditions[2],
# 	mar = c(2, 2, 4, 4), number.cex = 0.5, 
# 	number.digits = 2, tl.col = 'black',
# 	tl.cex = 1)
# dev.off()

# #######################
save.image(file = "healthy_pulp_20191206.Rdata")

