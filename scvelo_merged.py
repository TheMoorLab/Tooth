# scvelo_merged

import scvelo as scv
import os.path

# -----------------------------------------------------------------------------------
# Merged Periodontium 1, 2, 3 and 6
sdata 						= scv.read('/Users/laura/Documents/Data/perio/bustools_spliced_unspliced/mergednamed/sf.loom', cache=True)
ndata 						= scv.read('/Users/laura/Documents/Data/perio/bustools_spliced_unspliced/mergednamed/uf.loom', cache=True)
# testdata 					= scv.read('/Users/laura/Documents/Data/Peri/spliced.loom', cache=True) # we use this to get the umap coords/clusters. 
sdata.layers['spliced'] 	= sdata.X
sdata.layers['unspliced'] 	= ndata.X
sdata

scv.utils.show_proportions(sdata)

# import umap embeddings from R (Seurat) and add to anndata object


data_folder 			= "/Users/laura/Documents/Data/perio/bustools_spliced_unspliced"
# import harmony umap embedding saved from Seurat object
obsm  					= scv.read_csv(os.path.join(data_folder, "embednonames.csv"))
sdata.obsm["X_umap"] 	= obsm.values
scv.pp.filter_and_normalize(sdata, min_shared_counts=30, n_top_genes=2000) # don't run if using seurat pre-filtered data

scv.pp.moments(sdata, n_pcs=30, n_neighbors=30)
# scv.tl.umap(sdata) # if re-doing umap within scvelo
scv.tl.velocity(sdata)
scv.tl.velocity_graph(sdata)
scv.pl.velocity_embedding_stream(sdata, basis='umap')
scv.pl.velocity_embedding_grid(sdata, basis='umap')