library(tidyverse)
library(Seurat)
library(harmony)
library(cowplot)

#setwd("~/Desktop/EAE-Carmen/")

x = readRDS("1-integrated.rds")
#x$new <- ifelse(x$orig.ident == "Cont", "Ctrl", "PKD")
#x$orig.ident = NULL


# Linear dimensional reduction (PCA):

x <- ScaleData(x)
x = x  %>% RunPCA(npcs = 30)

#print(x[["pca"]], dims = 1:4, nfeatures = 5)
# first 4 principal components; for each PC, find top 5 genes contributing to the PC

#VizDimLoadings(x, dims = 1:2, reduction = "pca")

DimPlot(x, reduction = "pca", pt.size = 0.01) + facet_wrap(~ident)

ElbowPlot(x) # Determine the ‘dimensionality’ of the dataset for UMAP, heatmap, etc

#DimHeatmap(x, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(x, dims = 1:15, cells = 500, balanced = TRUE)


## UMAP before homony

#x <- x %>% FindNeighbors(reduction = "pca") %>% FindClusters(resolution = 0.5)# default resolution 0.8
#table(Idents(x))

x <- RunUMAP(x, dims = 1:10, verbose = F)

DimPlot(x, reduction = "umap", label = TRUE, pt.size = 0.01)
DimPlot(x, reduction = "umap", pt.size = 0.01) + facet_wrap(~ident) + NoLegend()


# run harmony

x <- x %>% RunHarmony("orig.ident")

#DimPlot(object = x, reduction = "harmony", pt.size = .4, group.by = "orig.ident")
#VlnPlot(object = x, features = "harmony_1", group.by = "orig.ident",  pt.size = .4)
#DimHeatmap(object = x, reduction = "harmony", cells = 500, dims = 1:3)

x <- x %>% RunUMAP(reduction = "harmony",  dims = 1:10)

saveRDS(x, file = "2-umap.rds", compress = TRUE)
x = readRDS("2-umap.rds")


# identify clusters

x <- x %>% FindNeighbors(reduction = "harmony") %>% FindClusters(resolution=0.35)
table(Idents(x))

DimPlot(x,reduction="umap", label=TRUE, pt.size = 0.01)# + NoLegend()

DimPlot(x,reduction="umap", group.by="orig.ident", pt.size = 0.01)

DimPlot(x,reduction="umap", label=TRUE, pt.size = 0.01, split.by = 'orig.ident', ncol = 2) 

FeaturePlot(x, "Dpp10", pt.size = 0.01, raster=FALSE)