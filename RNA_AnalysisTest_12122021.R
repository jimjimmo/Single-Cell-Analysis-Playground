"Chapters
1. Package Imports
2. Data Imports
3. Data QC and Inspection
5. Data Normalization
6. Data Clustering (PCA/UMAP)
7. Markers Identification
8. Putting all together
"
setwd("C:/Users/jimbo/Documents/R/Single-Cell-Analysis-Playground")

#Load Packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "C:/Users/jimbo/Documents/R")
#Initialize the Seurat object with the raw (non-normalized data)

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc

#Add columns using pbmc[[]]. Double bracket is useful for adding information
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

p1 <- FeatureScatter(pbmc, feature1 = "nFeature_RNA", feature2 = "nCount_RNA")
p2 <- FeatureScatter(pbmc, feature1 = "nFeature_RNA", feature2 = "percent.mt")
p3 <-FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
p1 + p2 + p3

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize")

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

#These tare the top 10 most variable genes
top10 <- head(VariableFeatures(pbmc), 10) #See the top 10 genes most interesting
top10 #display gene names

#Plot features with and without labels
p4 <- VariableFeaturePlot(pbmc)
p5 <- LabelPoints(plot = p4, points = top10, repel = TRUE)
p5

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

print(pbmc[["pca"]], dims = 1:2, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, nfeatures = 15, reduction = "pca") #showing the PC's and important genes

DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

#pbmc <- JackStraw(pbmc, num.replicate = 100)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

#JackStrawPlot(pbmc, dims = 1:15)

ElbowPlot(pbmc) #Figure out how many PCs needed to make clusters. In this case ~10 is best
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

#Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

pbmc <- RunUMAP(pbmc, dims = 1:20)
DimPlot(pbmc, reduction = "umap")