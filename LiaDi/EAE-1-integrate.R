library(tidyverse)
library(Seurat)
#setwd("~/Desktop/EAE-Carmen/")

# read from data folder
cre0 <- Read10X(data.dir = "CRE0/")
cre0 <- CreateSeuratObject(counts = cre0, 
                           project = "CRE0", 
                           assay = "RNA", 
                           min.cells = 100, 
                           min.features = 200)
head(cre0@meta.data)
table(cre0@meta.data$orig.ident)

cre1 <- Read10X(data.dir = "CRE1/")
cre1 <- CreateSeuratObject(counts = cre1, 
                           project = "CRE1", 
                           assay = "RNA", 
                           min.cells = 100, 
                           min.features = 200)
head(cre1@meta.data)

ncar0 <- Read10X(data.dir = "NCAR0/")
ncar0 <- CreateSeuratObject(counts = ncar0, 
                           project = "NCAR0", 
                           assay = "RNA", 
                           min.cells = 100, 
                           min.features = 200)
head(ncar0@meta.data)

ncar1 <- Read10X(data.dir = "NCAR1/")
ncar1 <- CreateSeuratObject(counts = ncar1,
                            project = "NCAR1", 
                            assay = "RNA", 
                            min.cells = 100, 
                            min.features = 200)
head(ncar1@meta.data)

x <- merge(ncar0, y=c(ncar1, cre0, cre1), add.cell.ids=c("NCAR0","NCAR1","CRE0","CRE1"), project="EAE")
table(x@meta.data$orig.ident)

# QC and selecting cells

rownames(x)[str_detect(rownames(x), "^mt")]

x[["percent.mt"]] = PercentageFeatureSet(x, pattern = "^mt-")

# 1)
VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 2) FeatureScatter: typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

p1 <- FeatureScatter(x, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1) #+ NoLegend()
p2 <- FeatureScatter(x, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.1) #+ NoLegend()
p1 + p2

# 3) filter
x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10)


# 4) Filter genes to match total features in the merged data (x)
# modified from Bing results; subset features present only in x
filter_genes <- function(seurat_obj) {
  seurat_obj <- subset(seurat_obj, features = rownames(x))
  return(seurat_obj)
}
x.list <- SplitObject(x)
x.list <- lapply(x.list, filter_genes)
rm(x)


# Integration

for (i in 1:length(x.list)) {
  x.list[[i]] <- NormalizeData(x.list[[i]], verbose = FALSE)
  x.list[[i]] <- FindVariableFeatures(x.list[[i]], selection.method="vst",nfeatures=2000,verbose =FALSE)
}
#rm(i)

# Warning: it takes really long: 1.5h
anchor <- FindIntegrationAnchors(object.list = x.list)
saveRDS(anchor, file = "1-anchor.rds", compress = TRUE)
rm(x.list)

x <- IntegrateData(anchorset = anchor)
rm(anchor)

saveRDS(x, file = "1-integrated.rds", compress = TRUE)
