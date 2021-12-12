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
pbmc.data <- Read10X(data.dir = "C:/Users/jimbo/Documents/R/Single-Cell-Analysis-Playground")
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

pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

cluster1.markers = FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25) #See which gene is a marker for cluster 1
head(cluster1.markers, n = 5) #See whichtop 5 genes that define this cluster

VlnPlot(pbmc, features = c(row.names(cluster1.markers)[1], row.names(cluster1.markers)[2]))

#Find all markers distinguishing cluster 5 from  clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0,3), min.pct = 0.25)
head(cluster5.markers, n = 5)

VlnPlot(pbmc, features = c(row.names(cluster5.markers)[1], row.names(cluster5.markers)[2])) #Shows that 5 and 1 may be very similar, but 5 and 3 and distinct

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

x <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
p7 <- FeaturePlot(pbmc, features = x$gene[1:4])
p8 <- FeaturePlot(pbmc, features = x$gene[5:8])


