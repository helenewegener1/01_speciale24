library(spaceRATScaffolds)
library(spaceRAT)
library(tidyverse)

################################################################################
################################## Load Data ################################### 
################################################################################

load(file='/home/people/s165827/data/Breast cancer classification/Bordet_paired_data/Bordet_RNA_gene_tpm.Rdata')
load(file="/home/people/s165827/data/Breast cancer classification/gtex.Rdata")
gtex_class_sample <- class_sample 
names(gtex_class_sample) <- colnames(gtex_scaffold)
gtex_scaffold_TPD_scaled <- readRDS('rds/gtex_scaffold_TPD_scaled.rds')

gtex_pheno <- as.data.frame(class_sample)
rownames(gtex_pheno) <- colnames(gtex_scaffold)

### Subset ###
gtex_groups <- readRDS(file="rds/gtex_groups.rds")
gtex_groups <- gtex_groups %>% unlist()
gtex_scaffold_subset <- gtex_scaffold[rownames(gtex_scaffold) %in% gtex_groups ,]

spaceRAT_gtex_scaffold_TPD_scaled_scaffold <- buildScaffold(gtex_scaffold_subset, 
                                                            pheno = gtex_pheno,
                                                            colname = 'class_sample',
                                                            data = 'exprs', 
                                                            subset_deg = FALSE,
                                                            threshold = 10,
                                                            # add_umap = TRUE,
                                                            pca_scale = TRUE)

plotScaffold(spaceRAT_gtex_scaffold_TPD_scaled_scaffold, dims = c(1,2), title = 'gtex scaffold top 20 DE genes') + theme(legend.position = "none")
# ggsave('plots/gtex_scaffold_DEsubset_dim12.png', width = 10, height = 8.72)
plotScaffold(spaceRAT_gtex_scaffold_TPD_scaled_scaffold, dims = c(2,3), title = 'gtex scaffold top 20 DE genes') + theme(legend.position = "none")
# ggsave('plots/gtex_scaffold_DEsubset_dim23.png', width = 10, height = 8.72)
plotScaffold(spaceRAT_gtex_scaffold_TPD_scaled_scaffold, dims = c(1,3), title = 'gtex scaffold top 20 DE genes') + theme(legend.position = "none")
# ggsave('plots/gtex_scaffold_DEsubset_dim13.png', width = 10, height = 8.72)


### Define centroids ###

# Simple spearman correlation to centroid (based on mean)
# centroids <- do.call(cbind, lapply(unique(gtex_class_sample), function(c) {rowMeans(gtex_scaffold_subset[,gtex_class_sample==c])}))
# colnames(centroids) <- unique(gtex_class_sample)
# cors <- do.call(rbind, lapply(1:584, function(i) {cor(gtex_scaffold_subset[,i], centroids, method = 'spearman')}))
# rownames(cors) <- colnames(gtex_scaffold_subset)

# Rank

## Define centroids 
centroids <- list()
for (i in unique(gtex_class_sample)) {
  centroids[[i]] <- rowMeans(gtex_scaffold_subset[,gtex_class_sample==i])
  
}

## Get distance to centroid for each sample
sample_distance_df <- c()
nc_pred <- c()
for(i in 1:ncol(gtex_scaffold_subset)) {
  centroid_distances <- c()
  for(ii in 1:length(centroids)) {
    centroid_distances <- c(centroid_distances, dist(rbind(gtex_scaffold_subset[,i], centroids[[ii]])))
    # centroid_distances <- c(centroid_distances, DistancePair(gtex_scaffold_subset[,i], centroids[[ii]]))
  }
  names(centroid_distances) <- names(centroids)
  # closest_samples_df[[colnames(gtex_scaffold_subset)[i]]] <- centroid_distances
  sample_distance_df <- rbind(sample_distance_df, c(colnames(gtex_scaffold_subset)[i], centroid_distances))
  nc_pred <- c(nc_pred, names(which.min(centroid_distances)))
}

sample_distance_df <- as.data.frame(sample_distance_df)
sample_distance_df <- column_to_rownames(sample_distance_df, 'V1')

## Define 10 samples closest to centroid 
closest_samples <- c()

for (tissue in unique(gtex_class_sample)){
  
  # tissue <- 'Blood'
  
  tissue_samples <- names(gtex_class_sample[which(gtex_class_sample == tissue)])
  
  distances <- sample_distance_df[rownames(sample_distance_df) %in% tissue_samples, tissue]
  sample_names <- rownames(sample_distance_df[rownames(sample_distance_df) %in% tissue_samples, ])
  distances <- as.double(distances)
  names(distances) <- sample_names
  
  top_10_closest_samples <- distances %>%
    sort() %>%
    head(n = 10) %>%
    names()
  
  closest_samples <- c(closest_samples, top_10_closest_samples)
  
}

### Subset for top 10 closest samples from full scaffold ###
gtex_scaffold_top10samples <- gtex_scaffold[, colnames(gtex_scaffold) %in% closest_samples]
gtex_pheno_top10samples <- gtex_pheno[rownames(gtex_pheno) %in% closest_samples ,]
names(gtex_pheno_top10samples) <- closest_samples
gtex_pheno_top10samples <- as.data.frame(gtex_pheno_top10samples)
colnames(gtex_pheno_top10samples) <- c('class_sample')

spaceRAT_gtex_scaffold_TPD_scaled_scaffold <- buildScaffold(gtex_scaffold_top10samples, 
                                                            pheno = gtex_pheno_top10samples,
                                                            colname = 'class_sample',
                                                            data = 'exprs', 
                                                            subset_deg = FALSE,
                                                            threshold = 10,
                                                            # add_umap = TRUE,
                                                            pca_scale = TRUE)

plotScaffold(spaceRAT_gtex_scaffold_TPD_scaled_scaffold, dims = c(1,2), title = 'gtex scaffold top 10 sample') + theme(legend.position = "none")
ggsave('plots/gtex_scaffold_top10sample_dim12.png', width = 10, height = 8.72)
plotScaffold(spaceRAT_gtex_scaffold_TPD_scaled_scaffold, dims = c(2,3), title = 'gtex scaffold top 10 sample') + theme(legend.position = "none")
ggsave('plots/gtex_scaffold_top10sample_dim23.png', width = 10, height = 8.72)
plotScaffold(spaceRAT_gtex_scaffold_TPD_scaled_scaffold, dims = c(1,3), title = 'gtex scaffold top 10 sample') + theme(legend.position = "none")
ggsave('plots/gtex_scaffold_top10sample_dim13.png', width = 10, height = 8.72)


### Subset for top 10 closest samples from subset (DE 20) scaffold ###
gtex_scaffold_subset_top10samples <- gtex_scaffold_subset[, colnames(gtex_scaffold_subset) %in% closest_samples]
gtex_pheno_top10samples <- gtex_pheno[rownames(gtex_pheno) %in% closest_samples ,]
names(gtex_pheno_top10samples) <- closest_samples
gtex_pheno_top10samples <- as.data.frame(gtex_pheno_top10samples)
colnames(gtex_pheno_top10samples) <- c('class_sample')

spaceRAT_gtex_scaffold_TPD_scaled_scaffold <- buildScaffold(gtex_scaffold_subset_top10samples, 
                                                            pheno = gtex_pheno_top10samples,
                                                            colname = 'class_sample',
                                                            data = 'exprs', 
                                                            subset_deg = FALSE,
                                                            threshold = 10,
                                                            # add_umap = TRUE,
                                                            pca_scale = TRUE)



plotScaffold(spaceRAT_gtex_scaffold_TPD_scaled_scaffold, dims = c(1,2), title = 'gtex scaffold top 20 DE genes + top 10 sample') + theme(legend.position = "none")
ggsave('plots/gtex_scaffold_DEsubset_top10sample_dim12.png', width = 10, height = 8.72)
plotScaffold(spaceRAT_gtex_scaffold_TPD_scaled_scaffold, dims = c(2,3), title = 'gtex scaffold top 20 DE genes + top 10 sample') + theme(legend.position = "none")
ggsave('plots/gtex_scaffold_DEsubset_top10sample_dim23.png', width = 10, height = 8.72)
plotScaffold(spaceRAT_gtex_scaffold_TPD_scaled_scaffold, dims = c(1,3), title = 'gtex scaffold top 20 DE genes + top 10 sample') + theme(legend.position = "none")
ggsave('plots/gtex_scaffold_DEsubset_top10sample_dim13.png', width = 10, height = 8.72)



################################################################################
################################# Project data #################################
################################################################################

# Bordet PCA
projectSample(scaffold = spaceRAT_gtex_scaffold_TPD_scaled_scaffold,
              sample = Bordet_RNA_gene_tpm, 
              title = 'PCA: Bordet samples projected onto GTEX scaffold.',
              dims = c(1, 3))

loadingPlot(spaceRAT_gtex_scaffold_TPD_scaled_scaffold, dims = c(1, 3))

ggsave('plots/Bordet_projected_on_gtex_genesubsetN_scaffold.png',
       width = 16,
       height = 10)

