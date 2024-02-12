library(Seurat)
library(tidyverse)
# library(clustree)

#######################  Load data ####################### 
seurat_object <- readRDS('rds/SingelCellBreastCancer_seurat_object.rds')

# Code from article 
# https://github.com/Swarbricklab-code/BrCa_cell_atlas/blob/main/scSubtype/Calculatingscoresandplotting.R

###################################################################
####################### Mitochondrial genes ####################### 
###################################################################

# Define mitochondrial genes 
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, 
                                                      pattern = "^MT.")
# Plot distribution 
seurat_object[[]] %>% 
  ggplot(aes(x = percent.mt)) +
  geom_density() + 
  geom_vline(xintercept = 20, color = 'red')

# Filter accoring to plot 
seurat_object <- subset(seurat_object, 
                        subset = percent.mt < 20)

###################################################################
######################### Ribosomal genes ######################### 
###################################################################

seurat_object <- seurat_object[!grepl('^RP[SL]', rownames(seurat_object)), ]

###################################################################
######################## Bad quality cells ######################## 
###################################################################

FeatureScatter(seurat_object, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') + geom_smooth(method = 'lm')
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 8000)

###################################################################
######################## Standard workflow ######################## 
###################################################################

seurat_object <- NormalizeData(object = seurat_object)
seurat_object <- FindVariableFeatures(object = seurat_object, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(x = seurat_object)
seurat_object <- ScaleData(object = seurat_object, features = all.genes)
seurat_object <- RunPCA(object = seurat_object, features = VariableFeatures(object = seurat_object))

###################################################################
########################### Clustering ############################ 
###################################################################

ElbowPlot(seurat_object)

seurat_object <- FindNeighbors(seurat_object, dims = 1:16)

seurat_object <- FindClusters(seurat_object, resolution = 0.5)
seurat_object <- FindClusters(seurat_object, resolution = 0.3)
seurat_object <- FindClusters(seurat_object, resolution = 0.2)
seurat_object <- FindClusters(seurat_object, resolution = 0.1)

seurat_object <- RunUMAP(seurat_object, dims = 1:16)

# clustree(seurat_object, prefix = "RNA_snn_res.")

Idents(seurat_object) <- "RNA_snn_res.0.2"

DimPlot(seurat_object, 
        reduction = "umap") + theme(legend.position = "none")

# Save object
saveRDS(seurat_object, 'rds/SingelCellBreastCancer_seurat_object_QC.rds')



