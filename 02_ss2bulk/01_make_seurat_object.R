library(Seurat)
library(glue)

data_dir <- '/home/people/s165827/data/singelCell_breast'
sample_dirs <- c('CID3586', 'CID3838', 'CID3921', 
                 'CID4040', 'CID4290A', 'CID4471',
                 'CID4465', 'CID44991', 'CID4513')

sample_subtypes <- c('HER2/ER_1', 'HER2_2', 'HER2_3',
                     'ER_1', 'ER_2', 'ER_3',
                     'TNBC_1', 'TNBC_2', 'TNBC_3')


# Make Seurat object per sample and save in list
seurat_object_list <- list()
for (sample in sample_dirs){
  
  expression_matrix <- ReadMtx(
    mtx = glue('{data_dir}/{sample}/count_matrix_sparse.mtx'), 
    features = glue('{data_dir}/{sample}/count_matrix_genes.tsv'),
    cells = glue('{data_dir}/{sample}/count_matrix_barcodes.tsv'),
    feature.column = 1,
  )
  
  seurat_object <- CreateSeuratObject(counts = expression_matrix,
                                      min.cells = 3)
  
  seurat_object_list[[sample]] <- seurat_object

}

# Combine list of Seurat objects into one Seurat object
seurat_object <- merge(seurat_object_list[[1]], 
                       y = seurat_object_list[2:length(seurat_object_list)], 
                       add.cell.ids = sample_subtypes, 
                       project = "SingelCellBreastCancer")

# seurat_object <- JoinLayers(seurat_object)

# head(seurat_object[[]])

# Save object
saveRDS(seurat_object, 'rds/SingelCellBreastCancer_seurat_object.rds')

# unique(sapply(X = strsplit(colnames(seurat_object), split = "_"), FUN = "[", 1))



