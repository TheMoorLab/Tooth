# Part 2 of calculating velocity
# Import loom files created by merging pulps 2,3,13 and 14 in Seurat
# Import umap embeddings from seurat to visualize velocity in the same umap coordinates

import scvelo as scv
import os.path

# -----------------------------------------------------------------------------------
# Open loom files from merged Pulp 2,3,13,14
# I ran python locally so make sure to check for correct paths
pulp_sdata 						= scv.read('/Users/delaura/Documents/Tooth/AllPulp/merged/sf.loom', cache=True)
pulp_ndata 						= scv.read('/Users/delaura/Documents/Tooth/AllPulp/merged/uf.loom', cache=True)
pulp_sdata.layers['spliced'] 	= pulp_sdata.X
pulp_sdata.layers['unspliced'] 	= pulp_ndata.X
pulp_sdata

scv.utils.show_proportions(pulp_sdata)

# import umap embeddings from R (Seurat) and add to anndata object
data_folder 			= "/Users/delaura/Documents/Tooth/AllPulp/merged"
# import harmony umap embedding saved from Seurat object
obsm  					= scv.read_csv(os.path.join(data_folder, "embednonames.csv"))
pulp_sdata.obsm["X_umap"] 	= obsm.values
scv.pp.filter_and_normalize(pulp_sdata, min_shared_counts=30, n_top_genes=2000) # don't run if using seurat pre-filtered data

scv.pp.moments(pulp_sdata, n_pcs=30, n_neighbors=30)
# scv.tl.umap(pulp_sdata) # if re-doing umap within scvelo
scv.tl.velocity(pulp_sdata)
scv.tl.velocity_graph(pulp_sdata)
# scv.pl.velocity_embedding_stream(pulp_sdata, basis='umap') # I find the stream visuals to look decieving 
scv.pl.velocity_embedding_grid(pulp_sdata, basis='umap', save= "perio_grid.pdf", dpi = 150)

