# Jaccard distance plotted tih force-directed graph layout
library(plyr) #join
library(proxy) # more distances (ejaccard)
library(qgraph) # force-directed graph layout

CTasgns 			= read.csv("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/CellTypesDefinitions.csv", header = FALSE)
names(CTasgns)		= c("cluster", 'metaCellType','cellType', 'colors')

Idents(object = merged_harmony) <- "groups_bysize"
av.ex_top 			= AverageExpression(merged_harmony,  add.ident = 'condition', features =  merged_harmony@assays$SCT@var.features)

av.ex.SCT2 			= data.frame(t(av.ex_top$SCT))

cluster_dist2   	= dist(av.ex.SCT2 , method = "ejaccard")
dist_mat2 		    = as.matrix(cluster_dist2)

# Set up node legend labels for qgraph
# Sub row names' underscore with space
rownames(av.ex.SCT2) = gsub("_", " ", rownames(av.ex.SCT2))
vertex_names 		 = rownames(av.ex.SCT2)

# node graph labels
nodes_numbers 		= 1:length(vertex_names)
row.names(dist_mat2) = nodes_numbers
colnames(dist_mat2) 	= nodes_numbers

# Set up groups for qgraph
conditions 				= c("perio", "pulp")
cc_collapsed_ids    	= lapply(conditions, function(a) grep(pattern = a, x = vertex_names)) 
names(cc_collapsed_ids) = conditions

################################ Varying vertex sizes by proportions #########################################
# Estimate proportions for each cell type relative to pulp and perio totals (since total population sizes are quite different between perio and pulp)
melted_celltype_condition2         = cbind.data.frame(cellType =merged_harmony@meta.data$groups_bysize, condition=merged_harmony@meta.data$condition)
melted_celltype_by_condition2      = as.data.frame.matrix(table(melted_celltype_condition2))
cellType_normalized_by_condition2  = as.data.frame(t(melted_celltype_by_condition2)/rowSums(t(melted_celltype_by_condition2)))

celltype_props_by_sample <- cellType_normalized_by_condition2 %>%
mutate(ind = factor(row_number())) %>%
gather(variable, value, -ind)

# Gathering and mutating cause the sample (perio/pulp) names to be replaced with (1,2), so we create a table with the key and merge back with celltype_props_by_sample
key_cat		 			 = data.frame(sample = c('perio', 'pulp'), ind = c(1,2))
merged_props 			 = merge(celltype_props_by_sample, key_cat)[-1]#
# We're only interested in the label column with combined labels of cell type (variable) and sample
celltype_props_by_sample = cbind.data.frame(prop= merged_props$value, CTsample = paste(merged_props$variable, merged_props$sample))

################################################################################################################
# Order proportions according to vertex order in distance matrix
v_names 				 = data.frame(CTsample=vertex_names)
ordered_props 			 = join(v_names, celltype_props_by_sample)
# Because the range of cell type proportions is very wide and the node size variation would too large to plot, we compress the values through a ln transform 
compress_sizes 		     = log(1000*ordered_props$prop) # the 1000 scaling is just so there are no values betwee 0 and 1 which would lead to a negative value after the log transform. 

# Force-directed graph layout of distance matrix
# Nodes are normalized by cell type proportion (relative to perio or pulp sample sizes)
dist_mi <- 1-dist_mat2  # 1 as qgraph takes similarity matrices as input
pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/PseudobulkCorrs", 
	"PerioPulpEjaccard_FDGL_normalized_node_sizes.pdf"), width=16, height=16)
# qgraph(dist_mi, layout='spring', vsize=3, groups = test, legend = TRUE)
qgraph(dist_mi, layout='spring', vsize=compress_sizes, color = c(rep('blue',length(cc_collapsed_ids[1])), rep('red', length(cc_collapsed_ids[2]))),  nodeNames = vertex_names, groups = cc_collapsed_ids, legend.mode="style2")
dev.off()

# Force-directed graph layout of distance matrix (node sizes are all equal)
dist_mi <- 1-dist_mat2  # 1 as qgraph takes similarity matrices as input
pdf(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/PseudobulkCorrs", 
	"PerioPulpEjaccard_FDGL_fixed_node_sizes.pdf"), width=16, height=16)
# qgraph(dist_mi, layout='spring', vsize=3, groups = test, legend = TRUE)
qgraph(dist_mi, layout='spring', vsize=3, nodeNames = vertex_names, groups = cc_collapsed_ids, legend.mode="style2")
dev.off()
