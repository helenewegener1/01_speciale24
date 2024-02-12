library(Seurat)
library(tidyverse)
library(SingleR)
library(celldex)
library(scran)

seurat_object <- readRDS('rds/SingelCellBreastCancer_seurat_object_batch_corrected.rds')


############################################################################
######################### Distribution of metrics ######################### 
############################################################################

metrics <-  c("nCount_RNA", "nFeature_RNA", "percent.mt")

FeaturePlot(seurat_object, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

############################################################################
################ Annotation using canonical lineage markers ################ 
############################################################################

cluster0_conserved_markers <- FindConservedMarkers(seurat_object,
                                                   ident.1 = 0,
                                                   grouping.var = "condition",
                                                   only.pos = TRUE,
                                                   min.pct = 0.25,  
                                                   min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)

############################################################################
################################# SingleR ################################## 
############################################################################

# SingleR 

ref <- HumanPrimaryCellAtlasData()

# subset from Becker GitHub
types_to_use <- c(
  "DC",
  "Epithelial_cells",
  "B_cell",
  "Neutrophils",
  "T_cells",
  "Monocyte",
  "Endothelial_cells",
  "Neurons",
  "Macrophage",
  "NK_cell",
  "BM",
  "Platelets",
  "Fibroblasts",
  "Astrocyte",
  "Myelocyte",
  "Pre-B_cell_CD34-",
  "Pro-B_cell_CD34+",
  "Pro-Myelocyte")

ref <- ref[,(colData(ref)$label.main %in% types_to_use)]

seurat_object.counts <- GetAssayData(seurat_object, 
                                     assay = 'RNA', 
                                     layer = 'counts')

pred <- SingleR(test = seurat_object.counts,
                ref = ref,
                labels = ref$label.main,
                de.method = 'wilcox')

pred_cluster <- SingleR(test = seurat_object.counts,
                        ref = ref,
                        labels = ref$label.main,
                        de.method = 'wilcox',
                        clusters = seurat_object[[]]$RNA_snn_res.0.2)

# seurat_object$cell_type_singleR <- pred$labels[match(rownames(seurat_object@meta.data), rownames(pred))]
seurat_object$cell_type_singleR <- pred$labels[match(seurat_object@meta.data$RNA_snn_res.0.2, rownames(pred_cluster))]

DimPlot(seurat_object, 
        reduction = 'umap',
        group.by = "cell_type_singleR")


############################################################################
################ Annotation using published gene signatures ################ 
############################################################################

gene_signatures <- c('EPCAM', 'MKI67', 'CD3C', 'CD68', 'MS4A1', 'JCHAIN', 'PECAM1', 'PDGFRB')
# epithelial cells (EPCAM)
# proliferating cells (MKI67)
# T cells (CD3D)
# myeloid cells (CD68)
# B cells (MS4A1)
# plasmablasts (JCHAIN)
# endothelial cells (PECAM1)
# mesenchymal cells (fibroblasts/perivascular-like cells; PDGFRB)

FeaturePlot(seurat_object, 
            reduction = "umap", 
            features = gene_signatures, 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE) 


temp_markers <-
  c("COL1A1","PDGFRA", "FAP", "PDPN", "CXCL12", "PDGFRB", "ACTA2",
    "CD36", "MCAM", "MYH11", "RGS5",
    "PECAM1", "VWF", "CD34", "ACKR1", "DLL4", "LYVE1", "MKI67")


FeaturePlot(seurat_object, 
            reduction = "umap", 
            features = temp_markers, 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE) 
