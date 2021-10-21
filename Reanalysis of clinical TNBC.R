library(Seurat)
library(harmony)
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


setwd("~/Analysis/STmeta")


data <- read.delim2("GSE118389_tpm_rsem.txt", header = TRUE)
TNBC <- CreateSeuratObject(counts = data, project = "TNBC")
TNBC <- NormalizeData(TNBC)
TNBC[["percent.mt"]] <- PercentageFeatureSet(TNBC, pattern = "^MT-")
TNBC <- FindVariableFeatures(TNBC, selection.method = "vst")
TNBC <- ScaleData(TNBC, verbose = FALSE)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ST.combined <- CellCycleScoring(TNBC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
TNBC <- RunPCA(TNBC, npcs = 100, ndims.print = 1:5)
DimPlot(TNBC, reduction = "pca")
DimHeatmap(TNBC, dims = 1:15, cells = 500, balanced = TRUE)


#######################
TNBC <- RunUMAP(TNBC, dims = 1:25, min.dist = 0.5)
TNBC <- RunTSNE(TNBC, dims = 1:15, check_duplicates = FALSE)
TNBC <- FindNeighbors(TNBC, dims = 1:15)
TNBC <- FindClusters(TNBC, resolution = 1)

DimPlot(TNBC, reduction = "tsne", pt.size = 2)

saveRDS(TNBC, file = "TNBC.rds")



ALL.markers <- FindAllMarkers(TNBC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "wilcox")
write.table(ALL.markers, file = "ALLmarkers.csv", sep = ",")


#######################
## Subset Epithelial cancer cell
EPI <- subset(TNBC, ident = c(2,3,5,7,8))
EPI <- RunUMAP(EPI, dims = 1:25, min.dist = 0.5)
EPI <- RunTSNE(EPI, dims = 1:15, check_duplicates = FALSE)
EPI <- FindNeighbors(EPI, dims = 1:15)
EPI <- FindClusters(EPI, resolution = 1)

DimPlot(EPI, pt.size = 2.5)
DimPlot(EPI, reduction = "tsne", pt.size = 2.5)

EPI.markers <- FindAllMarkers(EPI, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "wilcox")
write.table(EPI.markers, file = "EPImarkers.csv", sep = ",")

# 0 CCL20
# 1 CDH2
# 2 TBX3
# 3 RNASE1
# 4 PDGFA
# 5 CCL28
saveRDS(EPI, file = "EPI.rds")


### Module analysis
# use the gene signatures
# list1 is HMGA1 signatures
# list2 is CD44 signatures
list <- read.csv("list1.csv", header = FALSE)
list = unlist(list)
allgene = rownames(EPI)
list1 = list(allgene[allgene %in% unlist(list)])

list2 <- read.csv("list2.csv", header = FALSE)
list2 = unlist(list2)
allgene = rownames(EPI)
list2 = list(allgene[allgene %in% unlist(list2)])


mod1 <- AddModuleScore(EPI, features = list1, pool = allgene, name = "HMGA1sig.")
mod1 <- AddModuleScore(mod1, features = list2, pool = allgene, name = "CD44sig.")

write.table(as.matrix(mod1@meta.data[,7:8]), file = "signature.csv", sep = ",")

FeaturePlot(mod1, features = "HMGA1sig.1", cols = c("grey", "orange","red2") , pt.size = 2)
FeaturePlot(mod1, features = "CD44sig.1", cols = c("grey", "orange","red2") , pt.size = 2)
VlnPlot(mod1, features = "HMGA1sig.1")
VlnPlot(mod1, features = "CD44sig.1")
VlnPlot(subset(mod1, ident = c(2,3,5,7,8)), features = "HMGA1sig.1", pt.size = 2)
FeatureScatter(mod1, feature1 = "HMGA1sig.1", feature2 = "CD44sig.1", pt.size = 2)
FeaturePlot(mod1, features = c("HMGA1sig.1", "CD44sig.1"), blend = TRUE, pt.size = 2.5)

