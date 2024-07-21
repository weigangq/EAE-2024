#snrna merge workflow for cre0 and cre1 
#cre knockout 
#0 before signs of disease
#1 after signs of disease

library(SingleCellExperiment)
library(tidyverse)
library(scPred)
library(Seurat)
library(magrittr)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(gridExtra)
library(SingleR)
library(Seurat)
library(celldex)
library(pheatmap)
library(EnhancedVolcano)
library(scDblFinder)
library('beepr')
set.seed(123)



#load in data cre0

cre0 <- readRDS("/CRE0_preprocessed.rds")
cre1 <- readRDS("/cre1_preprocessed.rds")

#Merge 
#merging cre0 and cre1 

cre_merged <- merge(x = cre0, y = cre1, 
                    project = "Integrated_1272_1315")

head(cre_merged@meta.data)


#create a metadta data frame

metadata_cre <- cre_merged@meta.data
metadata_cre$cells <- rownames(metadata_cre)#add cell ids to metadata_cre_1315a
metadata_cre$sample <- NA
metadata_cre$sample[which(str_detect(metadata_cre$cells, "^cre0"))] <- "CRE0"
metadata_cre$sample[which(str_detect(metadata_cre$cells, "^cre1"))] <- "CRE1"

#rename columns 
metadata_cre <- metadata_cre %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# add meta data back to seurat obj

cre_merged@meta.data <- metadata_cre



metadata_cre %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")


#saving merged

saveRDS(cre_merged, file = "/Users/michellefranco/scRNA_patientss/preprocessed/cre_merged.rds")

#normalize merged

cre_integrated<- cre_merged
cre_integrated <- NormalizeData(cre_integrated)
cre_integrated <- FindVariableFeatures(cre_integrated)
cre_integrated <- ScaleData(cre_integrated)
cre_integrated <- RunPCA(cre_integrated)
cre_integrated <- FindNeighbors(cre_integrated)
cre_integrated <- FindClusters(cre_integrated, resolution = 0.3, cluster.name = "unintegrated_clusters")
ElbowPlot(cre_integrated)


#PCA Results

#DimPlot(cre_integrated, reduction = "pca") + NoLegend()
DimPlot(cre_integrated,
        reduction = "pca",
        label = TRUE,
        label.size = 6)

#UMAP

cre_integrated <- RunUMAP(cre_integrated, dims = 1:15, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(cre_integrated, reduction = "umap.unintegrated", group.by = c("sample", "unintegrated_clusters"))


####Intergration

# Plot the UMAP

DimPlot(cre_integrated,
        #reduction = "umap",
        label = TRUE,
        label.size = 6)

# UMAP of cells in each cluster by sample


DimPlot(cre_integrated, 
        label = TRUE, 
        split.by = "sample")  + NoLegend()


#extract identity and sample info

n_cells_cre0 <- FetchData(cre_integrated, 
                             vars = c("ident", "sample")) %>%
  dplyr::count(ident, sample)
n_cells_cre0

# Barplot of number of cells per cluster by sample
ggplot(n_cells_cre0, aes(x=ident, y=n, fill=sample)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_text(aes(label=n), vjust = -.2, position=position_dodge(1))


##############################################
#scpred
reference <- scPred::cre0

reference <- reference %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)

reference <- getFeatureSpace(reference, "cell_type")

reference <- trainModel(reference)

get_probabilities(reference) %>% head()
get_scpred(reference)
plot_probabilities(reference)








#comparing to unsupervised clustering


tab_cre0 <- table(Assigned=results_ptcre0$labels, Clusters=cre_integrated$seurat_clusters)
tab_cre0


pheatmap(log10(tab_cre0+10), color = colorRampPalette(c('white','blue'))(10))


#DEG between samples 

Idents(cre_integrated) <- cre_integrated@meta.data$sample


deResults_cre0 <- FindMarkers(cre_integrated,
                                 ident.1 = 'Normal', 
                                 ident.2 = 'Tumor', 
                                 min.pct = 0.25, 
                                 logfc.threshold = 0.25)
head(deResults_cre0)


ggplot(deResults_cre0, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.5) + 
  theme_minimal() + 
  labs(x = "avg_log2FC", y = "-Log10 p_val_adj", title = "Volcano Plot") +
  theme(legend.position = "e") 



deResults_cre0$GeneName <- rownames(deResults_cre0)



EnhancedVolcano(deResults_cre0,
                lab = deResults_cre0$GeneName,
                x = 'avg_log2FC',
                y = 'p_val_adj')




sessionInfo()


