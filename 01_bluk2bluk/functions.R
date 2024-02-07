library(glue)
library(patchwork)

subset_sample <- function(scaffold, class_sample, n_samples, rank = FALSE){
  
  # scaffold <- gtex_scaffold
  # n_samples <- 10
  # # rank <- TRUE
  
  names(class_sample) <- colnames(scaffold)
  
  if (rank == TRUE) {
    scaffold <- apply(scaffold, 2, rank)
  }
  
  centroids <- list()
  for (i in unique(class_sample)) {
    centroids[[i]] <- rowMeans(scaffold[,class_sample==i])
    
  }
  
  ## Get distance to centroid for each sample
  sample_distance_df <- c()
  nc_pred <- c()
  for(i in 1:ncol(scaffold)) {
    centroid_distances <- c()
    for(ii in 1:length(centroids)) {
      centroid_distances <- c(centroid_distances, dist(rbind(scaffold[,i], centroids[[ii]])))
      # if (rank == FALSE) {
      #   centroid_distances <- c(centroid_distances, dist(rbind(scaffold[,i], centroids[[ii]])))
      # } else if (rank == TRUE){
      #   centroid_distances <- c(centroid_distances, DistancePair(scaffold[,i], centroids[[ii]]))
      # }
    }
    names(centroid_distances) <- names(centroids)
    sample_distance_df <- rbind(sample_distance_df, c(colnames(scaffold)[i], centroid_distances))
    nc_pred <- c(nc_pred, names(which.min(centroid_distances)))
  }
  
  sample_distance_df <- as.data.frame(sample_distance_df)
  sample_distance_df <- column_to_rownames(sample_distance_df, 'V1')
  
  ## Define 10 samples closest to centroid 
  closest_samples <- c()
  
  for (tissue in unique(class_sample)){
    
    # tissue <- 'Blood'
    
    tissue_samples <- names(class_sample[which(class_sample == tissue)])
    
    distances <- sample_distance_df[rownames(sample_distance_df) %in% tissue_samples, tissue]
    sample_names <- rownames(sample_distance_df[rownames(sample_distance_df) %in% tissue_samples, ])
    distances <- as.double(distances)
    names(distances) <- sample_names
    
    top_closest_samples <- distances %>%
      sort() %>%
      head(n = n_samples) %>%
      names()
    
    closest_samples <- c(closest_samples, top_closest_samples)
    
  }
  
  # Make scaffold data
  scaffold_top_samples <- scaffold[, colnames(scaffold) %in% closest_samples]

  # Make pheno 
  pheno <- as.data.frame(class_sample)
  pheno_top_samples <- pheno[rownames(pheno) %in% closest_samples ,] %>% as.data.frame()
  rownames(pheno_top_samples) <- colnames(scaffold_top_samples)
  colnames(pheno_top_samples) <- c('class_sample')
  
  # Return 
  return_list <- list(scaffold_top_samples, pheno_top_samples)
  
  return(return_list)
  
  
}

subset_scaffold <- function(basline_scaffold, class_sample, rds_path, extra_name = '') {
  
  # Read rds file 
  DEgenes <- readRDS(rds_path)
  name <- gsub("rds/|\\.rds", "", rds_path)
  
  # Scaffild
  DEgenes <- DEgenes %>% unlist()
  scaffold_DEgenes_subset <- basline_scaffold[rownames(basline_scaffold) %in% DEgenes ,]
  
  # Pheno
  pheno <- as.data.frame(class_sample)
  rownames(pheno) <- colnames(scaffold_DEgenes_subset)
  
  DEgenes_Scaffold <- buildScaffold(scaffold_DEgenes_subset, 
                                    pheno = pheno,
                                    colname = 'class_sample',
                                    data = 'exprs', 
                                    subset_deg = TRUE,
                                    threshold = 10,
                                    # add_umap = TRUE,
                                    pca_scale = TRUE)
  
  
  p1 <- plotScaffold(DEgenes_Scaffold, dims = c(1,2), title = glue('DE genes ({name}) PC1 & PC2 {extra_name}')) + theme(legend.position = "none")
  ggsave(glue('plots/gtex_scaffold_{name}_dim12{extra_name}.png'), width = 10, height = 8.72)
  p2 <- plotScaffold(DEgenes_Scaffold, dims = c(2,3), title = glue('DE genes ({name}) PC2 & PC3 {extra_name}')) + theme(legend.position = "none")
  ggsave(glue('plots/gtex_scaffold_{name}_dim23{extra_name}.png'), width = 10, height = 8.72)
  p3 <- plotScaffold(DEgenes_Scaffold, dims = c(1,3), title = glue('DE genes ({name}) PC1 & PC3 {extra_name}')) + theme(legend.position = "none")
  ggsave(glue('plots/gtex_scaffold_{name}_dim13{extra_name}.png'), width = 10, height = 8.72)
  
  return_plot <- wrap_plots(p1, p2, p3, ncol = 2)
  ggsave(glue('plots/gtex_scaffold_{name}_wrapped{extra_name}.png'), width = 15, height = 8)
  
  return(return_plot)
  
}


