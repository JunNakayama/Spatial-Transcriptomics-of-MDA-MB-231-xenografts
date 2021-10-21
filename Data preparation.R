library(Seurat)
library(ComplexHeatmap)
library(scater)
library(Matrix)
library(tidyverse)
library(dplyr)
library(tidyr)
library(reshape2)
library(pheatmap)
library(Rtsne)
library(cowplot)
library(ggplot2)
library(ggsci)
library(scales)
library(MAST)
library(DOSE)
library(patchwork)
library(plotly)
library(monocle)
library(MASS)
library(loomR)
library(RColorBrewer)
library(grDevices)
library(colorRamps)
library(data.table)
library(hexbin)




## For human cancer cell analysis

PT3hTPM <- read.csv("PT3hTPM.csv",header=TRUE)
hDF = data.frame(PT3hTPM)
hpDF = subset(hDF, gene_biotype == c("protein_coding"))

hpDFdata <- data.matrix(hpDF[,-c(1:4)])

PTs = hpDFdata[,1:93]
LYs = hpDFdata[,94:137]


#Seurat Object preparation
rownames(PTs) = paste(hpDF$gene_name)
PTs = data.frame(PTs)
rownames(LYs) = paste(hpDF$gene_name)
LYs = data.frame(LYs)

STPTs <- CreateSeuratObject(PTs, project = "SPATIAL_TRANSCRIPTOMICS")
STLYs <- CreateSeuratObject(LYs, project = "SPATIAL_TRANSCRIPTOMICS")

STPTs@meta.data$stim <- "PT"
STLYs@meta.data$stim <- "Lymph"
ST.combined <- merge(STPTs, y = c(STLYs))


#######################
## Filtering 

ST.combined <- NormalizeData(ST.combined)
ST.combined[["percent.mt"]] <- PercentageFeatureSet(ST.combined, pattern = "MT.")

ST.combined <- FindVariableFeatures(ST.combined, selection.method = "vst")
ST.combined <- ScaleData(ST.combined, verbose = FALSE, vars.to.regress = "percent.mt")


## Cell Cycle scoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ST.combined <- CellCycleScoring(ST.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

ST.combined <- RunPCA(ST.combined, npcs = 100, ndims.print = 1:5)
DimPlot(ST.combined, reduction = "pca")
DimHeatmap(ST.combined, dims = 1:15, cells = 500, balanced = TRUE)
VizDimLoadings(ST.combined, dims = 1:2, reduction = "pca")


ST.combined <- JackStraw(ST.combined, num.replicate = 100)
ST.combined <- ScoreJackStraw(ST.combined, dims = 1:20)

JackStrawPlot(ST.combined, dims = 1:15)
ElbowPlot(ST.combined)


#######################
###　Clustering and UMAP plot
ST.combined <- FindNeighbors(ST.combined, dims = 1:15)
ST.combined <- FindClusters(ST.combined, resolution = 1)
ST.combined <- RunUMAP(ST.combined, dims = 1:25, min.dist = 0.5)
ST.combined <- RunTSNE(ST.combined, dims = 1:15, check_duplicates = FALSE)

DimPlot(ST.combined, reduction = "tsne", pt.size = 2, group.by = "stim")
saveRDS(ST.combined, file = "ST.combined.rds")

#######################
### find marker genes
ST.com.markers <- FindAllMarkers(ST.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "wilcox")
write.table(ST.com.markers, file = "STallmarkers.csv", sep = ",")




#######################
#######################
## For mouse stromal cell analysis

PT3mTPM <- read.csv("PT3mTPM.csv",header=TRUE)

mDF = data.frame(PT3mTPM)
mpDF = subset(mDF, gene_biotype == c("protein_coding"))

mpDFdata <- data.matrix(mpDF[,-c(1:4)])

mPTs = mpDFdata[,1:93]
mLYs = mpDFdata[,94:137]



#Seurat Object preparation
rownames(mPTs) = paste(mpDF$gene_name)
mPTs = data.frame(mPTs)
rownames(mLYs) = paste(mpDF$gene_name)
mLYs = data.frame(mLYs)

mSTPTs <- CreateSeuratObject(mPTs, project = "SPATIAL_TRANSCRIPTOMICS", min.cells = 5)
mSTLYs <- CreateSeuratObject(mLYs, project = "SPATIAL_TRANSCRIPTOMICS", min.cells = 5)

mSTPTs@meta.data$stim <- "PT"
mSTLYs@meta.data$stim <- "Lymph"
mST.combined <- merge(mSTPTs, y = c(mSTLYs))

#######################
## Filtering 

mST.combined <- NormalizeData(mST.combined)
mST.combined[["percent.mt"]] <- PercentageFeatureSet(mST.combined, pattern = "Mt.")

mST.combined <- FindVariableFeatures(mST.combined, selection.method = "vst")
mST.combined <- ScaleData(mST.combined, verbose = FALSE, vars.to.regress = "percent.mt")


## Cell Cycle scoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
mST.combined <- CellCycleScoring(mST.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


mST.combined <- RunPCA(mST.combined, npcs = 100, ndims.print = 1:5)
DimPlot(mST.combined, reduction = "pca")
DimHeatmap(mST.combined, dims = 1:15, cells = 500, balanced = TRUE)
VizDimLoadings(mST.combined, dims = 1:2, reduction = "pca")


mST.combined <- JackStraw(mST.combined, num.replicate = 100)
mST.combined <- ScoreJackStraw(mST.combined, dims = 1:20)

JackStrawPlot(mST.combined, dims = 1:15)
ElbowPlot(mST.combined)


#######################
###　Clustering and UMAP plot
mST.combined <- FindNeighbors(mST.combined, dims = 1:15)
mST.combined <- FindClusters(mST.combined, resolution = 1)
mST.combined <- RunUMAP(mST.combined, dims = 1:25, min.dist = 0.5)
mST.combined <- RunTSNE(mST.combined, dims = 1:15, check_duplicates = FALSE)

DimPlot(mST.combined, reduction = "tsne", pt.size = 2, group.by = "stim")
saveRDS(mST.combined, file = "mST.combined.rds")


#######################
### find marker genes
mST.com.markers <- FindAllMarkers(mST.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "wilcox")
write.table(mST.com.markers, file = "mSTallmarkers.csv", sep = ",")


