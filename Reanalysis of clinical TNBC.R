library(Seurat)
library(harmony)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(MAST)



setwd("~/Analysis/STmeta")


####### GSE163210


T1 <-  Read10X_h5("~/Analysis/STmeta/Revise/GSE163210_RAW/GSM4974849_T1.h5")
T2 <-  Read10X_h5("~/Analysis/STmeta/Revise/GSE163210_RAW/GSM4974849_T1.h5")
T3 <-  Read10X_h5("~/Analysis/STmeta/Revise/GSE163210_RAW/GSM4974851_T3.h5")
CTC1 <-  Read10X_h5("~/Analysis/STmeta/Revise/GSE163210_RAW/GSM4974852_CTC1.h5")
CTC2 <-  Read10X_h5("~/Analysis/STmeta/Revise/GSE163210_RAW/GSM4974853_CTC2.h5")

T1 <- CreateSeuratObject(T1, project = "GSE163210_T1")
T2 <- CreateSeuratObject(T2, project = "GSE163210_T2")
T3 <- CreateSeuratObject(T3, project = "GSE163210_T3")
CTC1 <- CreateSeuratObject(CTC1, project = "GSE163210_CTC1")
CTC2 <- CreateSeuratObject(CTC2, project = "GSE163210_CTC2")


GSE163210 <- merge(T1, y = c(T2, T3, CTC1, CTC2))


GSE163210 <- NormalizeData(GSE163210, normalization.method = "LogNormalize", scale.factor = 10000)
GSE163210 <- FindVariableFeatures(GSE163210)
GSE163210[["percent.mt"]] <- PercentageFeatureSet(GSE163210, pattern = "MT-")
GSE163210 <- ScaleData(GSE163210, vars.to.regress = "percent.mt")
GSE163210 <- RunPCA(GSE163210, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
GSE163210 <- FindNeighbors(GSE163210, reduction = "pca", dims = 1:75, nn.eps = 0.5)
GSE163210 <- FindClusters(GSE163210, resolution = 3, n.start = 10)
GSE163210 <- RunUMAP(GSE163210, dims = 1:75, min.dist = 0.75)

DimPlot(GSE163210)

saveRDS(GSE163210, file = "GSE163210.rds")



FeaturePlot(GSE163210, features = c("HMGA1", "CD44"), blend = TRUE)




# > GSE163210
# An object of class Seurat 
# 54328 features across 8494 samples within 1 assay 
# Active assay: RNA (54328 features, 2000 variable features)
#  2 dimensional reductions calculated: pca, umap


gene <- c("HMGA1", "CD44")
obj <- GSE163210[unlist(gene),]
mat <- t(as.matrix(GetAssayData(obj)))
write.table(mat, file = "HMGA1CD44-PDX.csv", sep = ",")


write.table(GSE163210@meta.data, file = "GSE163210meta.csv", sep = ",")









#### Integrated dataset
GSE118389 <- readRDS("./GSE118389.rds")
GSE161529 <- readRDS("./GSE161529.rds")
GSE176078 <- readRDS("./GSE176078.rds")
GSE180286 <- readRDS("./GSE180286.rds")
Yamamotolab <- readRDS("./Yamamotolab.rds")

D <- merge(GSE118389, y = c(GSE161529, GSE176078, GSE180286, Yamamotolab), merge.data = TRUE)


metadata <- read.csv("meta-predata.csv")
metadata <- data.frame(metadata)
rownames(metadata) <- metadata[,1]

D@meta.data <- metadata

cl = c("Primary", "Metastasis", "DCIS")
D2 <- subset(D, subset = Tissue_state %in% cl)

D2 <- NormalizeData(D2)
VlnPlot(D2, features = "percent.mt", group.by = "Dataeset")
D2 <- subset(D2, subset = nFeature_RNA > 500 & percent.mt < 20)


D2 <- FindVariableFeatures(D2)
D2 <- ScaleData(D2, vars.to.regress = "percent.mt")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
D2 <- CellCycleScoring(D2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)



D2 <- RunPCA(D2, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
D2 <- RunHarmony(D2, "orig.ident")

D2 <- RunUMAP(D2, reduction = "harmony", dims = 1:100, min.dist = 0.5)
D2 <- FindNeighbors(D2, dims = 1:15, reduction = "harmony")
D2 <- FindClusters(D2, resolution = 1)



saveRDS(D2, file = "BRCAonly.rds")




### select the TNBC patients from three datasets

TNBC <- subset(BRCAonly, subset = subtype == "TNBC")

TNBC <- RunUMAP(TNBC, reduction = "harmony", dims = 1:100, min.dist = 0.5)
TNBC <- FindNeighbors(TNBC, dims = 1:15, reduction = "harmony")
TNBC <- FindClusters(TNBC, resolution = 1)

DimPlot(TNBC)

saveRDS(TNBC, file = "TNBC.rds")






sel = c(0, 3, 6, 7, 11, 12, 15, 16, 19, 20, 21, 23, 24)
EPI = subset(TNBC, ident = sel)

EPI <- RunUMAP(EPI, reduction = "harmony", dims = 1:100, min.dist = 0.5)
EPI <- FindNeighbors(EPI, dims = 1:15, reduction = "harmony")
EPI <- FindClusters(EPI, resolution = 1)


DimPlot(EPI)
DimPlot(EPI, group.by = "Tissue_state")




FeaturePlot(EPI, "HMGA1")
FeaturePlot(EPI, "CD44")
FeaturePlot(EPI, features = c("HMGA1", "CD44"), blend = TRUE)
FeatureScatter(EPI, feature1 = "HMGA1", feature2 = "CD44", group.by = "orig.ident")








####### GSE163210

