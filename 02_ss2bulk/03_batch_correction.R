library(Seurat)
library(tidyverse)

# Read rds
seurat_object <- readRDS('rds/SingelCellBreastCancer_seurat_object_QC.rds')

# Add column with subtype
pattern <- "^[^_]*"
seurat_object[["condition"]] <- regmatches(
  rownames(seurat_object@meta.data), 
  regexpr(pattern, rownames(seurat_object@meta.data))
)

# Add column with sample id
pattern <- "^[^_]*_[^_]*"
seurat_object[["sample_id"]] <- regmatches(
  rownames(seurat_object@meta.data), 
  regexpr(pattern, rownames(seurat_object@meta.data))
)

# See batch effect
DimPlot(seurat_object, 
        reduction = "umap",
        split.by = "condition")

########################################################################
####################### Correct for batch effect ####################### 
########################################################################

options(future.globals.maxSize = 1000 * 1024^2)

seurat_object <- IntegrateLayers(object = seurat_object, 
                                 method = RPCAIntegration,
                                 orig.reduction = "pca", 
                                 new.reduction = 'integrated.rpca',
                                 reference = 1, 
                                 verbose = FALSE)

seurat_object[["RNA"]] <- JoinLayers(seurat_object[["RNA"]])

########################################################################
############################ Refind clusters ########################### 
########################################################################

DefaultAssay(object = seurat_object) <- "RNA"

seurat_object <- ScaleData(seurat_object, verbose = FALSE)
seurat_object <- RunPCA(seurat_object, npcs = 50, verbose = FALSE)
seurat_object <- FindNeighbors(object = seurat_object, dims = 1:30)
seurat_object <- FindClusters(object = seurat_object, resolution = 0.2)
seurat_object <- RunUMAP(seurat_object, reduction = "pca", dims = 1:30)

# Check batch effect now
DimPlot(seurat_object,
        reduction = "umap",
        split.by = "condition")

Idents(seurat_object) <- "RNA_snn_res.0.2"

# Save object
saveRDS(seurat_object, 'rds/SingelCellBreastCancer_seurat_object_batch_corrected.rds')

