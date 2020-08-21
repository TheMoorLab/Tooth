library(Seurat) # 3.1.3
library(ggplot2)
library(sctransform)
library(harmony)

#prior code seurat V2: 20190523_perio_updated.R


# perio1 <- Read10X(data.dir = "/IMCR_shares/Moorlab/Common/Tooth_project/CellRanger_output/perio/v2_C92_scRNAseq_4_Perio1")
# perio1 <- CreateSeuratObject(counts = perio1,   project = "perio1",min.cells = 3,min.features = 200)

# perio2 <- Read10X(data.dir = "/IMCR_shares/Moorlab/Common/Tooth_project/CellRanger_output/perio/v2_C92_scRNAseq_4_Perio2")
# perio2 <- CreateSeuratObject(counts = perio2,   project = "perio2",min.cells = 3,min.features = 200)

# perio3 <- Read10X(data.dir = "/IMCR_shares/Moorlab/Common/Tooth_project/CellRanger_output/perio/v2_scRNAseq3_Period_WisToothAdultMildCaries")
# perio3 <- CreateSeuratObject(counts = perio3,   project = "perio3",min.cells = 3,min.features = 200)

# perio4 <- Read10X(data.dir = "/IMCR_shares/Moorlab/Common/Tooth_project/CellRanger_output/perio/v2_scRNAseq_6_Carious_Periodontium")
# perio4 <- CreateSeuratObject(counts = perio6,   project = "perio6",min.cells = 3,min.features = 200)

# perio5 <- Read10X(data.dir = "/IMCR_shares/Moorlab/Common/Tooth_project/CellRanger_output/perio/v2_4_scRNAseq_WisdomTooth_Upper_NoTotErupt_7_Perio")
# perio5 <- CreateSeuratObject(counts = perio7,   project = "perio7",min.cells = 3,min.features = 200)


perio1 <- Read10X(data.dir = "/IMCR_shares/Moorlab/Common/Tooth_project/CellRanger_output/perio/Perio1")
perio1 <- CreateSeuratObject(counts = perio1,   project = "perio1",min.cells = 3,min.features = 200)

perio2 <- Read10X(data.dir = "/IMCR_shares/Moorlab/Common/Tooth_project/CellRanger_output/perio/Perio2")
perio2 <- CreateSeuratObject(counts = perio2,   project = "perio2",min.cells = 3,min.features = 200)

perio3 <- Read10X(data.dir = "/IMCR_shares/Moorlab/Common/Tooth_project/CellRanger_output/perio/Perio3")
perio3 <- CreateSeuratObject(counts = perio3,   project = "perio3",min.cells = 3,min.features = 200)

perio4 <- Read10X(data.dir = "/IMCR_shares/Moorlab/Common/Tooth_project/CellRanger_output/perio/Perio4")
perio4 <- CreateSeuratObject(counts = perio4,   project = "perio4",min.cells = 3,min.features = 200)

perio5 <- Read10X(data.dir = "/IMCR_shares/Moorlab/Common/Tooth_project/CellRanger_output/perio/Perio5")
perio5 <- CreateSeuratObject(counts = perio5,   project = "perio5",min.cells = 3,min.features = 200)

all_perio 		<- merge(x=perio1, y=c(perio2,perio3,perio4,perio5),add.cell.ids=c("perio1","perio2","perio3","perio4","perio5"))

table(all_perio$orig.ident)

# store mitochondrial percentage in object meta data
all_perio 		<- PercentageFeatureSet(all_perio, pattern = "^MT-", col.name = "percent.mt")

# QC measures prior to subsetting
VlnPlot(all_perio, features = c("percent.mt"), 
        pt.size = 0.2, ncol = 4) + theme(axis.title = element_text(size = 20), axis.text= element_text(size = 20)) 
ggsave(filename = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPerio/PostReview/", "mt_presubset.pdf"),  width=14, height=10)

VlnPlot(all_perio, features = c("nCount_RNA"), 
        pt.size = 0.2, ncol = 4) + theme(axis.title = element_text(size = 20), axis.text= element_text(size = 20)) 
ggsave(filename = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPerio/PostReview/", "nCount_presubset.pdf"),  width=14, height=10)

VlnPlot(all_perio, features = c("nFeature_RNA"), 
        pt.size = 0.2, ncol = 4) + theme(axis.title = element_text(size = 20), axis.text= element_text(size = 20)) 
ggsave(filename = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPerio/PostReview/", "nFeature_presubset.pdf"),  width=14, height=10)

# Backup Seurat object and save

all_perio_backup <- all_perio
all_perio 		 <- subset(all_perio, subset = percent.mt < 20 & nCount_RNA < 50000)

# QC measures AFTER subsetting
VlnPlot(all_perio, features = c("percent.mt"), 
        pt.size = 0.2, ncol = 4) + theme(axis.title = element_text(size = 20), axis.text= element_text(size = 20)) 
ggsave(filename = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPerio/PostReview/", "mt_postsubset.pdf"),  width=14, height=10)

VlnPlot(all_perio, features = c("nCount_RNA"), 
        pt.size = 0.2, ncol = 4) + theme(axis.title = element_text(size = 20), axis.text= element_text(size = 20)) 
ggsave(filename = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPerio/PostReview/", "nCount_postsubset.pdf"),  width=14, height=10)

VlnPlot(all_perio, features = c("nFeature_RNA"), 
        pt.size = 0.2, ncol = 4) + theme(axis.title = element_text(size = 20), axis.text= element_text(size = 20)) 
ggsave(filename = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPerio/PostReview/", "nFeature_postsubset.pdf"),  width=14, height=10)


# # run sctransform
# all_perio <- SCTransform(all_perio, vars.to.regress = "percent.mt", verbose = TRUE)

# save(all_perio,file = "/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/r_output/all_perio_20191204.rds")

# #seurat pipeline
# all_perio <- RunPCA(all_perio, verbose = FALSE)
# all_perio <- RunUMAP(all_perio, dims = 1:50, verbose = FALSE)
# all_perio <- FindNeighbors(all_perio, dims = 1:50, verbose = FALSE)
# all_perio <- FindClusters(all_perio, verbose = FALSE)
# DimPlot(all_perio, label = TRUE) + NoLegend()

# all_perio$clustering <- Idents(object = all_perio)
# Idents(object = all_perio) <- "clustering"
# Idents(object = all_perio) <- "orig.ident"

# #harmony
# set.seed(20191204)
# all_perio<-RunHarmony(all_perio,"orig.ident", theta = 2, plot_convergence = TRUE, nclust = 50, 
#                         max.iter.cluster = 20, max.iter.harmony = 5)
# set.seed(20191204)
# all_perio <- RunUMAP(all_perio, dims = 1:50, verbose = FALSE,reduction = "harmony")
# all_perio <- FindNeighbors(all_perio, dims = 1:50, verbose = FALSE,reduction = "harmony")
# all_perio <- FindClusters(all_perio, verbose = FALSE,reduction="harmony")

# DimPlot(all_perio, label = TRUE,reduction = "umap")

# FeaturePlot(all_perio,features=c("THY1","MYH11","COL1A1","KRT14","PECAM1","PTPRC","MBP","GFRA3","SPOCK3","MS4A1"), reduction = "umap",order = T)
# FeaturePlot(all_perio,features=c("MYH11"), reduction = "umap",order = T)

# all_markers<-FindAllMarkers(all_perio)

# save(all_perio,all_markers,file = "/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/r_output/all_perio_markers_20191204.rds")
# save.image(file = "perio_20191204.Rdata")




