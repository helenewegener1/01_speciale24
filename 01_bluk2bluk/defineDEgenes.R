

load(file="/home/people/s165827/data/Breast cancer classification/gtex.Rdata")
gtex_pheno <- as.data.frame(class_sample)
rownames(gtex_pheno) <- colnames(gtex_scaffold)

library(genefilter) # You may need to load this library for the coef function

log2fold_changes <- data.frame(row.names = rownames(gtex_scaffold))

labels <- class_sample

# Find top DE genes
DE20groups <- list()
DE_L2FC_0.3_pval_0.05 <- list()
DE_L2FC_0.2_pval_0.05 <- list()
DE_L2FC_0.1_pval_0.05 <- list()
DE_pval_0.05 <- list()

for(i in unique(labels)) {
  
  pvals <- apply(gtex_scaffold, 1, function(x) wilcox.test(x[labels == i], x[!labels == i])$p.value)
  pvals <- p.adjust(pvals)
  
  # Calculate log2fold changes and filter
  for (j in 1:ncol(gtex_scaffold)) {
    group1 <- gtex_scaffold[labels == i, j]  # Expression values for group i
    group2 <- gtex_scaffold[labels != i, j]  # Expression values for other groups

    # Calculate log2fold change using the mean expression values
    log2fold_change <- log2(mean(group1, na.rm = TRUE)) - log2(mean(group2, na.rm = TRUE))

    # Store the log2fold change in the dataframe
    log2fold_changes[j, i] <- log2fold_change

  }
  
  # Append significant features based on adjusted p-values and log2fold changes
  significant_indices_L2FC_0.3_pval_0.05 <- which(abs(log2fold_changes[, i]) > 0.3 & pvals < 0.05)
  significant_indices_L2FC_0.2_pval_0.05 <- which(abs(log2fold_changes[, i]) > 0.2 & pvals < 0.05)
  significant_indices_L2FC_0.1_pval_0.05 <- which(abs(log2fold_changes[, i]) > 0.1 & pvals < 0.05)
  significant_indices_pval_0.05 <- which(pvals < 0.05)


  DE_L2FC_0.3_pval_0.05[[i]] <- c(na.omit(names(pvals)[significant_indices_L2FC_0.3_pval_0.05]))
  DE_L2FC_0.2_pval_0.05[[i]] <- c(na.omit(names(pvals)[significant_indices_L2FC_0.2_pval_0.05]))
  DE_L2FC_0.1_pval_0.05[[i]] <- c(na.omit(names(pvals)[significant_indices_L2FC_0.1_pval_0.05]))
  DE_pval_0.05[[i]] <- c(na.omit(names(pvals)[significant_indices_pval_0.05]))
  
  top_20_significant_indices <- sort(pvals)[1:20]
  DE20groups[[i]] <- names(top_20_significant_indices)
  
}

saveRDS(DE_L2FC_0.3_pval_0.05, 'rds/DE_L2FC_0.3_pval_0.05.rds')
saveRDS(DE_L2FC_0.2_pval_0.05, 'rds/DE_L2FC_0.2_pval_0.05.rds')
saveRDS(DE_L2FC_0.1_pval_0.05, 'rds/DE_L2FC_0.1_pval_0.05.rds')
saveRDS(DE_pval_0.05, 'rds/DE_pval_0.05.rds')
saveRDS(DE20groups, 'rds/20DEgenes.rds')


# # De normaliserede ekspressionsværdier udregnes således:
# gtex_scaffold_DE_scaled <- NULL
# for(i in 1:length(groups)) {
#   pca <- prcomp(t(gtex_scaffold[groups[[i]],]), scale. = TRUE)
#   gtex_scaffold_DE_scaled <- rbind(gtex_scaffold_DE_scaled, as.matrix(gtex_scaffold[groups[[i]],]/pca$sdev[1]^2))
# }
# gtex_scaffold_DE_scaled
# 
# saveRDS(gtex_scaffold_DE_scaled, 'rds/gtex_scaffold_20DEgenes.rds')
# 
# gtex_scaffold_DE_scaled_rank <- apply(gtex_scaffold_DE_scaled, 2, rank)
# saveRDS(gtex_scaffold_DE_scaled_rank, 'rds/gtex_scaffold_20DEgenes_rank.rds')


############ TCGA ############ 
### find de genes ###

library(genefilter) # You may need to load this library for the coef function

log2fold_changes <- data.frame(row.names = rownames(tcga_scaffold))

labels <- class_sample

# Find top DE genes
DE20groups <- list()
DE_L2FC_0.3_pval_0.05 <- list()
DE_pval_0.05 <- list()

for(i in unique(labels)) {
  
  pvals <- apply(tcga_scaffold, 1, function(x) wilcox.test(x[labels == i], x[!labels == i])$p.value)
  pvals <- p.adjust(pvals)
  
  # Calculate log2fold changes and filter
  for (j in 1:ncol(tcga_scaffold)) {
    group1 <- tcga_scaffold[labels == i, j]  # Expression values for group i
    group2 <- tcga_scaffold[labels != i, j]  # Expression values for other groups
    
    # Calculate log2fold change using the mean expression values
    log2fold_change <- log2(mean(group1, na.rm = TRUE)) - log2(mean(group2, na.rm = TRUE))
    
    # Store the log2fold change in the dataframe
    log2fold_changes[j, i] <- log2fold_change
    
  }
  
  # Append significant features based on adjusted p-values and log2fold changes
  significant_indices_L2FC_0.3_pval_0.05 <- which(abs(log2fold_changes[, i]) > 0.3 & pvals < 0.05)
  significant_indices_pval_0.05 <- which(pvals < 0.05)
  
  DE_L2FC_0.3_pval_0.05[[i]] <- c(na.omit(names(pvals)[significant_indices_L2FC_0.3_pval_0.05]))
  DE_pval_0.05[[i]] <- c(na.omit(names(pvals)[significant_indices_pval_0.05]))
  
  top_20_significant_indices <- sort(pvals)[1:20]
  DE20groups[[i]] <- names(top_20_significant_indices)
  
}

saveRDS(DE_L2FC_0.3_pval_0.05, 'rds/tcga_DE_L2FC_0.3_pval_0.05.rds')
saveRDS(DE_pval_0.05, 'rds/tcga_DE_pval_0.05.rds')
saveRDS(DE20groups, 'rds/tcga_20DEgenes.rds')
