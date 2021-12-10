"Chapters
1. Package Imports
2. Data Imports
3. Data QC and Inspection
5. Data Normalization
6. Data Clustering (PCA/UMAP)
7. Markers Identification
8. Putting all together
"
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
