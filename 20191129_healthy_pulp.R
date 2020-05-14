library(Seurat)
library(ggplot2)
library(sctransform)
library(harmony)

#prior code seurat V2: 20190210_tooth_integrated_harmony_2.R


pulp2 <- Read10X(data.dir = "/IMCR_shares/Moorlab/Common/Tooth_project/CellRanger_output/pt2_adult_healthy/pulp")
pulp2 <- CreateSeuratObject(counts  = pulp2,   project = "pulp2",min.features = 200,min.cells = 3)

pulp3 <- Read10X(data.dir = "/IMCR_shares/Moorlab/Common/Tooth_project/CellRanger_output/sample3_adult_healthy/pulp")
pulp3 <- CreateSeuratObject(counts = pulp3,   project = "pulp3",min.features = 200,min.cells = 3)

pulp4 <- Read10X(data.dir = "/IMCR_shares/Moorlab/Common/Tooth_project/CellRanger_output/sample4_adult_healthy/pulp")
pulp4 <- CreateSeuratObject(counts = pulp4,   project = "pulp4",min.cells = 3,min.features = 200)

pulp13 <- Read10X(data.dir = "/IMCR_shares/Moorlab/Common/Tooth_project/CellRanger_output/scRNAseq_13_1_AdultPulpHealthy_14_24")
pulp13 <- CreateSeuratObject(counts = pulp13,   project = "pulp13v3",min.cells = 3,min.features = 200)

pulp14 <- Read10X(data.dir = "/IMCR_shares/Moorlab/Common/Tooth_project/CellRanger_output/scRNAseq_14_2_AdultPulpHealthy_2_35_45")
pulp14 <- CreateSeuratObject(counts = pulp14,   project = "pulp14v3",min.cells = 3,min.features = 200)

all_healthy<-merge(x=pulp2, y=c(pulp3,pulp4,pulp13,pulp14),add.cell.ids=c("pulp2","pulp3","pulp4","pulp13","pulp14"))

table(all_healthy$orig.ident)

# store mitochondrial percentage in object meta data
all_healthy <- PercentageFeatureSet(all_healthy, pattern = "^MT-", col.name = "percent.mt")


VlnPlot(all_healthy, features = c("percent.mt"), 
        pt.size = 0.2, ncol = 4)
VlnPlot(all_healthy, features = c("nCount_RNA"), 
        pt.size = 0.2, ncol = 4)

all_healthy_backup<-all_healthy

all_healthy<-subset(all_healthy, subset = percent.mt < 20 & nCount_RNA < 25000)


# run sctransform
all_healthy <- SCTransform(all_healthy, vars.to.regress = "percent.mt", verbose = TRUE)

save(all_healthy,file = "/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/r_output/all_healthy_20191129.rds")

#seurat pipeline
all_healthy <- RunPCA(all_healthy, verbose = FALSE)
all_healthy <- RunUMAP(all_healthy, dims = 1:50, verbose = FALSE)
all_healthy <- FindNeighbors(all_healthy, dims = 1:50, verbose = FALSE)
all_healthy <- FindClusters(all_healthy, verbose = FALSE)
DimPlot(all_healthy, label = TRUE) + NoLegend()

all_healthy$clustering <- Idents(object = all_healthy)
Idents(object = all_healthy) <- "clustering"
Idents(object = all_healthy) <- "orig.ident"

#harmony
set.seed(20191129)
all_healthy<-RunHarmony(all_healthy,"orig.ident", theta = 2, plot_convergence = TRUE, nclust = 50, 
                        max.iter.cluster = 20, max.iter.harmony = 5)
set.seed(20191129)
all_healthy <- RunUMAP(all_healthy, dims = 1:50, verbose = FALSE,reduction = "harmony")
all_healthy <- FindNeighbors(all_healthy, dims = 1:50, verbose = FALSE,reduction = "harmony")
all_healthy <- FindClusters(all_healthy, verbose = FALSE,reduction="harmony")

DimPlot(all_healthy, label = TRUE,reduction = "umap")

FeaturePlot(all_healthy,features=c("THY1","MYH11","COL1A1","KRT14","PECAM1","PTPRC","MBP","GFRA3","SPOCK3","MS4A1"), reduction = "umap",order = T)
FeaturePlot(all_healthy,features=c("STEAP4"), reduction = "umap",order = T)

all_markers<-FindAllMarkers(all_healthy)

save(all_healthy,all_markers,file = "/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/r_output/all_healthy_markers_20191129.rds")
save.image(file = "healthy_pulp_20191129.Rdata")

