# Updated 11/02/2020
# run scvelo of harmony-merged cell subset
library(data.table)
library(Matrix)
load("~/Documents/Data/perio/harmony_allperio_no15.rds") 

# IMPORTANT this seurat object "test" is from seurat2 and in seurat3 the S4 class has been renamed from seurat to Seurat. 
# This will cause issues when using methods from the old seurat therefore the object needs to be updated as below:
# seu = UpdateSeuratObject(test)
# for more info: https://github.com/satijalab/seurat/issues/990

# harmony_bcodes  = names(test@ident) 
harmony_bcodes  = merged_harmony@assays$RNA@data@Dimnames[[2]]
harmony_ids 	= numeric()
perios    		= c(1,2,3,6)
## perios    		= c(1,2,3,6,7)  # non matching barcodes in sample 7 so we're leaving it out for now...

# Local mac
# data_folder 	= "/Users/laura/Documents/Data/perio/bustools_spliced_unspliced"

# Science Cloud
data_folder 	 = "/IMCR_shares/Moorlab/Laura/MPro/Data/Perios/bustools_spliced_unspliced"

merged_spliced 	 = list()

merged_unspliced = list()

# for (p in 1:length(perios)) {
for (p in 2:length(perios)) {


	i 						= perios[p]
	# Find all of the current perio's bcodes 
	label 					= paste0('perio', i)
	# barcodes for current sample
	curr_bcodes 			= harmony_bcodes[harmony_bcodes %like% label]
	# ids in harmony (important in case we don't use all samples and therefore cells that went into harmony)
	harmony_ids 			= c(harmony_ids, which(harmony_bcodes %like% label))
	# trimmed barcodes to remove the perio prefix for matching to bcodes from bustools' outputs "spliced" and "unspliced"
	bcodes 					= gsub(".*_","",curr_bcodes)

	# SPLICED
	# Import list of all bcodes for that sample's spliced file
	label_all_s_bcodes    		= read.csv(file.path(data_folder, label,'spliced','s.barcodes.txt'), header = FALSE)
	names(label_all_s_bcodes) 	= 'spliced'
	# find id of all bcode matches in this perio's s.barcodes.txt files to find ids to use on s.mtx matrix
	curr_s_bids 				= which(label_all_s_bcodes$spliced %in% bcodes) 
	# write.csv(curr_s_bids, file = file.path(data_folder, label, 'spliced', paste0(label, "_bids.csv")), row.names = FALSE)
	# temp3						= Read10X(data.dir = file.path(data_folder, label,'spliced'))
	spliced_mat 				= readMM(file.path(data_folder, label, 'spliced', 's.mtx'))
	merged_spliced[[p]] 		= spliced_mat[curr_s_bids,] 

	# UNSPLICED
	# Same as above but let's find the cells in the unspliced file
	label_all_u_bcodes    		= read.csv(file.path(data_folder, label,'unspliced','u.barcodes.txt'), header = FALSE)
	names(label_all_u_bcodes) 	= 'unspliced'
	curr_u_bids					= which(label_all_u_bcodes$unspliced %in% bcodes) 
	# write.csv(curr_u_bids, file = file.path(data_folder, label, 'unspliced', paste0(label, "_bids.csv")), row.names = FALSE)

	unspliced_mat 				= readMM(file.path(data_folder, label, 'unspliced', 'u.mtx'))
	merged_unspliced[[p]] 		= unspliced_mat[curr_u_bids,] 
	# merged_bcodes[[p]]			= curr_bcodes

}

# merge un/spliced matrices across samples
mergeds 				= do.call(rbind, merged_spliced)
mergedu 				= do.call(rbind, merged_unspliced)

# orig_bcodes 			= names(test@ident[1:2520]) # that's to exclude last sample (perio7)
# orig_bcodes 			= harmony_bcodes[1:2520] # 

trimmed_bcodes 			= gsub(".*_","",orig_bcodes)

# mergedu@Dimnames[[2]] 	= names(tempb) 
# allgns 					= read.csv('/Users/laura/Documents/Data/perio/bustools_spliced_unspliced/perio1/spliced/s.genes.txt', header= FALSE)
# mergedu@Dimnames[[1]] 	= as.character(allgns$V1) 

# write merged barcodes
write.table(trimmed_bcodes, file=file.path(data_folder, 'merged', 'spliced', 's.barcodes.txt'), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
write.table(trimmed_bcodes, file=file.path(data_folder, 'merged', 'unspliced', 'u.barcodes.txt'), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)

# save merged spliced/unspliced matrices
writeMM(mergeds, file=file.path(data_folder, 'merged', 'spliced','s.mtx')) 
writeMM(mergedu, file=file.path(data_folder, 'merged', 'unspliced', 'u.mtx'))

# Check that writing of MM was sucessful
testS = readMM(file=file.path(data_folder, 'merged', 'spliced','s.mtx')) 
testU = readMM(file=file.path(data_folder, 'merged', 'unspliced', 'u.mtx'))

# ---------------------
# write merged barcodes including unique prefix per sample
write.table(orig_bcodes, file=file.path(data_folder, 'mergednamed', 'spliced', 's.barcodes.txt'), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
write.table(orig_bcodes, file=file.path(data_folder, 'mergednamed', 'unspliced', 'u.barcodes.txt'), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)

####################################################################################
# NOW go into Rstudio in server to use bustools and do the following (making sure all paths point to the right places)
# 
# c(spliced, unspliced) %<-% read_velocity_output(spliced_dir = "/mnt/bustools_velocity/mergednamed/spliced",
#                                                 spliced_name = "s",
#                                                 unspliced_dir = "/mnt/bustools_velocity/mergednamed/unspliced",
#                                                 unspliced_name = "u")
# 
# seu <- CreateSeuratObject(splicedn, assay = "sf")
# seu[["uf"]] <- CreateAssayObject(unspliced)
# 
## Choose server directory to save results in and save Seurat object
# outputsDir = "/mnt/bustools_velocity/mergednamed"
# save(seu, file = file.path(outputsDir,"seu_merged.RData")
####################################################################################


# Copy Seurat object to local dir
# BACK to local Rstudio
load('/Users/laura/Documents/Data/perio/bustools_spliced_unspliced/mergednamed/seu_merged.RData')

# srobj$RNA@meta.features needs to have columns otherwise loomR throws an error. Populate it running:
seu2 <- FindVariableFeatures(object = seu)

as.loom(seu2, assay = 'uf', filename = file.path(data_folder, 'mergednamed',"uf.loom"))
as.loom(seu2, assay = 'sf', filename = file.path(data_folder, 'mergednamed',"sf.loom"))

# Check you have the correct number of barcodes matching harmony numbers
tot_bcodes = sum(unlist(lapply(har_bids, length)))

seuHarmony = UpdateSeuratObject(test) # necessary for Seurat 2 objects which is the case of : ~/Documents/Data/perio/harmony_allperio_no15.rds

# Extract umap embeddings on ids we chose 
embeds = Embeddings(seuHarmony[["umap"]])[harmony_ids,]
# save umap embedding
write.csv(embeds, file = file.path(data_folder, "embednonames.csv"), row.names = FALSE)
