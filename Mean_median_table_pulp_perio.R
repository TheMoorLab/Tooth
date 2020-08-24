
outputsDir 					 = "/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses"
# PERIO
# Genes stats
genes_per_cell_perio 		 = Matrix::colSums(all_perio@assays$RNA@counts > 0) 
perio_genes_stat 			 = rbind(`median genes` = median(genes_per_cell_perio), `mean genes`= mean(genes_per_cell_perio), `sd genes` = sd(genes_per_cell_perio))
colnames(perio_genes_stat)   = 'Periodontium'

# Number of cells stats
perio_cellnum_stat 			 = rbind(`median number of cells` = median(table(all_perio@meta.data$orig.ident)), `mean number of cell`= mean(table(all_perio@meta.data$orig.ident)), `sd number of cells` = sd(table(all_perio@meta.data$orig.ident)))
colnames(perio_cellnum_stat) = 'Periodontium'

perio_stats = rbind(perio_genes_stat, perio_cellnum_stat )
# PULP
# Genes stats
genes_per_cell_pulp 		 = Matrix::colSums(all_healthy@assays$RNA@counts > 0) 
pulp_genes_stat 			 = rbind(`median genes` = median(genes_per_cell_pulp ), `mean genes`= mean(genes_per_cell_pulp ), `sd genes` = sd(genes_per_cell_pulp ))
colnames(pulp_genes_stat) 	 = 'Pulp'

# Number of cells stats
pulp_cellnum_stat 			 = rbind(`median number of cells` = median(table(all_healthy@meta.data$orig.ident)), `mean number of cell`= mean(table(all_healthy@meta.data$orig.ident)), `sd number of cells` = sd(table(all_healthy@meta.data$orig.ident)))
colnames(pulp_cellnum_stat)  = 'Pulp'

pulp_stats 					 = rbind(pulp_genes_stat, pulp_cellnum_stat )


write.table(cbind(perio_stats, pulp_stats) , file.path(outputsDir,'Perio_pulp_stats.csv'),
                                    col.names=TRUE, row.names=TRUE, sep=", ")

