# Updated 11/02/2020
# Prepare data to run scvelo of harmony-merged cell subset
# I had to use science cloud R instance http://172.23.177.146:8787 to use bustools 
library(data.table)
library(Matrix)

load("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/r_output/all_perio_markers_20191204.rds")

harmony_bcodes  = all_perio@assays$RNA@data@Dimnames[[2]]
harmony_ids 	= numeric()
perios    		= c(1,2,3,6)
## perios    		= c(1,2,3,6,7)  # non matching barcodes in sample 7 so we're leaving it out for now...

# Science Cloud
data_folder 	 = "/IMCR_shares/Moorlab/Laura/MPro/Data/Perios/bustools_spliced_unspliced"
merged_spliced 	 = list()
merged_unspliced = list()

for (p in 1:length(perios)) {

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

orig_bcodes 			= harmony_bcodes[1:2497] # that's to exclude last sample (perio7); total is otherwise 2883
trimmed_bcodes 			= gsub(".*_","",orig_bcodes)

output_dir <- file.path(file.path(data_folder, 'merged'))
if (!dir.exists(output_dir)){
dir.create(output_dir)
} else {
    print("Dir already exists!")
}

output_dir_spliced <- file.path(file.path(output_dir, 'spliced'))
if (!dir.exists(output_dir_spliced)){
dir.create(output_dir_spliced)
} else {
    print("Dir already exists!")
}

output_dir_unspliced <- file.path(file.path(output_dir, 'unspliced'))
if (!dir.exists(output_dir_unspliced)){
dir.create(output_dir_unspliced)
} else {
    print("Dir already exists!")
}
# write merged barcodes
write.table(trimmed_bcodes, file=file.path(output_dir_spliced, 's.barcodes.txt'), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
write.table(trimmed_bcodes, file=file.path(output_dir_unspliced, 'u.barcodes.txt'), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)

# save merged spliced/unspliced matrices
writeMM(mergeds, file=file.path(output_dir_spliced,'s.mtx')) 
writeMM(mergedu, file=file.path(output_dir_unspliced, 'u.mtx'))

# Check that writing of MM was sucessful
testS = readMM(file=file.path(output_dir_spliced,'s.mtx')) 
testU = readMM(file=file.path(output_dir_unspliced, 'u.mtx'))

# ---------------------
# write merged barcodes including unique prefix per sample
write.table(orig_bcodes, file=file.path(output_dir_spliced, 's.barcodes.txt'), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
write.table(orig_bcodes, file=file.path(output_dir_unspliced, 'u.barcodes.txt'), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)

####################################################################################
library(BUSpaRse) # to read bustools velocity output
library(zeallot) # For %<-% that unpacks lists in the Python manner

c(spliced, unspliced) %<-% read_velocity_output(spliced_dir = output_dir_spliced,
                                                spliced_name = "s",
                                                unspliced_dir = output_dir_unspliced,
                                                unspliced_name = "u")

seu <- CreateSeuratObject(spliced, assay = "sf")
seu[["uf"]] <- CreateAssayObject(unspliced)

# srobj$RNA@meta.features needs to have columns otherwise loomR throws an error. Populate it running:
seu2 <- FindVariableFeatures(object = seu)

as.loom(seu2, assay = 'uf', filename = file.path(output_dir, "uf.loom"))
as.loom(seu2, assay = 'sf', filename = file.path(output_dir, "sf.loom"))

# Extract umap embeddings on ids we chose 
embeds = Embeddings(all_perio[["umap"]])[harmony_ids,]

# save umap embedding
write.csv(embeds, file = file.path(output_dir, "embednonames.csv"), row.names = FALSE)

# save workspace
save.image(file = file.path(data_folder, "perio_vel.RData"))
