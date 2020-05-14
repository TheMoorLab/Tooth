

library(Seurat)
library(ggplot2)
library(sctransform)
library(harmony)

#prior code seurat V2: 20191111_perio_pulp_all.R


perio1 <- Read10X(data.dir = "/IMCR_shares/Moorlab/Common/Tooth_project/CellRanger_output/perio/v2_C92_scRNAseq_4_Perio1")
perio1 <- CreateSeuratObject(counts = perio1,   project = "perio1",min.cells = 3,min.features = 200)

perio2 <- Read10X(data.dir = "/IMCR_shares/Moorlab/Common/Tooth_project/CellRanger_output/perio/v2_C92_scRNAseq_4_Perio2")
perio2 <- CreateSeuratObject(counts = perio2,   project = "perio2",min.cells = 3,min.features = 200)

perio3 <- Read10X(data.dir = "/IMCR_shares/Moorlab/Common/Tooth_project/CellRanger_output/perio/v2_scRNAseq3_Period_WisToothAdultMildCaries")
perio3 <- CreateSeuratObject(counts = perio3,   project = "perio3",min.cells = 3,min.features = 200)

perio6 <- Read10X(data.dir = "/IMCR_shares/Moorlab/Common/Tooth_project/CellRanger_output/perio/v2_scRNAseq_6_Carious_Periodontium")
perio6 <- CreateSeuratObject(counts = perio6,   project = "perio6",min.cells = 3,min.features = 200)

perio7 <- Read10X(data.dir = "/IMCR_shares/Moorlab/Common/Tooth_project/CellRanger_output/perio/v2_4_scRNAseq_WisdomTooth_Upper_NoTotErupt_7_Perio")
perio7 <- CreateSeuratObject(counts = perio7,   project = "perio7",min.cells = 3,min.features = 200)


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




merged<-merge(x=perio1, y=c(perio2,perio3,perio6,perio7,pulp2,pulp3,pulp4,pulp13,pulp14),add.cell.ids=c("perio1","perio2","perio3","perio6","perio7","pulp2","pulp3","pulp4","pulp13","pulp14"))

table(merged$orig.ident)

# store mitochondrial percentage in object meta data
merged <- PercentageFeatureSet(merged, pattern = "^MT-", col.name = "percent.mt")


VlnPlot(merged, features = c("percent.mt"), 
        pt.size = 0.2, ncol = 4)
VlnPlot(merged, features = c("nCount_RNA"), 
        pt.size = 0.2, ncol = 4)

merged_backup<-merged

merged<-subset(merged, subset = percent.mt < 20 & nCount_RNA < 25000)


# run sctransform
merged <- SCTransform(merged, vars.to.regress = "percent.mt", verbose = TRUE)

save(merged,file = "/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/r_output/merged_20191204.rds")

#seurat pipeline
set.seed(20191204)
merged <- RunPCA(merged, verbose = FALSE)
set.seed(20191204)
merged <- RunUMAP(merged, dims = 1:50, verbose = FALSE)
merged <- FindNeighbors(merged, dims = 1:50, verbose = FALSE)
merged <- FindClusters(merged, verbose = FALSE)
DimPlot(merged, label = TRUE) + NoLegend()

merged$clustering <- Idents(object = merged)
Idents(object = merged) <- "clustering"
Idents(object = merged) <- "orig.ident"

#harmony
set.seed(20191204)
merged_harmony<-RunHarmony(merged,"orig.ident", theta = 2, plot_convergence = TRUE, nclust = 50, 
                      max.iter.cluster = 20, max.iter.harmony = 5)
set.seed(20191204)
merged_harmony <- RunUMAP(merged_harmony, dims = 1:50, verbose = FALSE,reduction = "harmony")
merged_harmony <- FindNeighbors(merged_harmony, dims = 1:50, verbose = FALSE,reduction = "harmony")
merged_harmony <- FindClusters(merged_harmony, verbose = FALSE,reduction="harmony")

DimPlot(merged_harmony, label = TRUE,reduction = "umap")
merged_harmony$clustering <- Idents(object = merged_harmony)
Idents(object = merged_harmony) <- "clustering"
Idents(object = merged_harmony) <- "orig.ident"

FeaturePlot(merged_harmony,features=c("THY1","MYH11","COL1A1","KRT14","PECAM1","PTPRC","MBP","GFRA3","SPOCK3","MS4A1"), reduction = "umap",order = T)
FeaturePlot(merged_harmony,features=c("percent.mt"), reduction = "umap",order = T)

all_markers<-FindAllMarkers(merged_harmony)

#currently markers not run yet
save(merged_harmony,file = "/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/r_output/merged_harmony_markers_20191204.rds")

save(merged_harmony,all_markers,file = "/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/r_output/merged_harmony_markers_20191204.rds")
save.image(file = "perio_20191204.Rdata")

