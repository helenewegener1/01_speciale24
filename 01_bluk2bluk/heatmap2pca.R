library(ggplot2)
library(mclust)
library(biomaRt)
library(patchwork)
library(genefilter)
library(hgu133plus2.db)
library(rankdist)
library(tidyverse)
library(cowplot)
library(pheatmap)
library(ComplexHeatmap)
# library(bioDist)

load("/home/people/s165827/data/Breast cancer classification/CITBCMST.rdata")
load("/home/people/s165827/data/Bordet_paired_data/Bordet.rdata")
# source('/home/people/s165827/data/Breast cancer classification/estimate_cbp.R')
# source('../hel_functions.R')

pam50_lvls <- c('Normal', 'LumA', 'LumB', 'Her2', 'Basal')
pam50_colors <- c("#00BA38", "#619CFF", "#00BFC4", "#B79F00", "#F8766D"); names(pam50_colors) <- pam50_lvls

cit_lvls <- c('normL', 'lumA', 'lumB', 'lumC', 'mApo', 'basL')
cit_colors <- c("#00BA38", "#619CFF", "#00BFC4", "#F564E3", "#B79F00", "#F8766D"); names(cit_colors) <- cit_lvls
sample_names <- c("HER2-03", "HER2-21", "LUMA-18", "LUMA-24", "LUMA-27", "LUMA-29", "LUMB-01", "LUMB-17", "TN-18", "TN-22"); names(sample_names) <- paste('sample', 1:10)

# # Subset core 
# Bordet_array_core <- Bordet_array[rownames(Bordet_array) %in% rownames(CIT_core), ]
# Bordet_RNA_core <- Bordet_RNA_tpm[rownames(Bordet_RNA_tpm) %in% rownames(CIT_core), ]
# 
# # Rank normalization 
# Bordet_array_rank <- apply(Bordet_array_core, 2, rank)
# Bordet_RNA_rank <- apply(Bordet_RNA_core, 2, rank)



################################################################################
################################# GTEx Heatmap ################################# 
################################################################################


# Load scaffold
load(file="/home/people/s165827/data/Breast cancer classification/gtex.Rdata")

testdata <- Bordet_RNA_tpm
rownames(testdata) <- Bordet_annot[match(rownames(testdata), Bordet_annot$Probe.Set.ID),]$Gene.Symbol
## Subset both sets based on genes in common

testdata <- testdata[intersect(rownames(testdata), rownames(gtex_scaffold)),]
# testdata <- testdata[, sample_names]
gtex_scaffold <- gtex_scaffold[intersect(rownames(testdata), rownames(gtex_scaffold)),]

## Rank both
gtex_scaffold_rank <- apply(gtex_scaffold,2,rank)
testdata_rank <- apply(testdata,2,rank)

# Simple spearman correlation to centroid (based on mean)
centroids <- do.call(cbind, lapply(unique(class_sample), function(c) {rowMeans(gtex_scaffold[,class_sample==c])}))
colnames(centroids) <- unique(class_sample)
cors <- do.call(rbind, lapply(1:10, function(i) {cor(testdata_rank[,i], centroids, method = 'spearman')}))
rownames(cors) <- names(sample_names)

# heatmap(cors)
# cors <- do.call(rbind, lapply(1:10, function(i) {cor(testdata[,i], gtex_scaffold, method = 'spearman')}))
# rownames(cors) <- paste('sample', 1:10)
library(ComplexHeatmap)
# column_ha = HeatmapAnnotation(GTEx = class_sample)
gtex_plot <- Heatmap(cors, name = "Spearman's\ncorrelation", show_column_names = T, column_title_rot = 90, cluster_rows = F)
#              cluster_column_slices = T, top_annotation = column_ha, column_split = class_sample)

# pdf('/Users/lro 1 2/Library/CloudStorage/GoogleDrive-larsronnolsen@gmail.com/My Drive/Research/Projects/Ongoing/Subtyping framework/Analysis/Breast cancer classification/Plots_CBP/Figure3_hm_gtex.pdf', height = 6, width = 8)
gtex_plot <- draw(gtex_plot, padding = unit(c(43.7, 2, 2, 2), "mm"), column_title = "GTEx heatmap")
# dev.off()
# 
pdf('../plots/Figure_3A.pdf', height = 8, width = 10)
draw(gtex_plot, padding = unit(c(43.7, 2, 2, 2), "mm"), column_title = "GTEx heatmap", column_title_side = 'top')
draw(gtex_plot, column_title = "GTEx heatmap")
dev.off()

################################################################################
################################# TCGA Heatmap ################################# 
################################################################################

# Load scaffold
load(file="/home/people/s165827/data/Breast cancer classification/tcga_scaffold_v2.Rdata")

# Import expression from samples (this will be in any type of format - pre-processing script will follow)
testdata <- Bordet_RNA_tpm
rownames(testdata) <- Bordet_annot[match(rownames(testdata), Bordet_annot$Probe.Set.ID),]$Gene.Symbol
# testdata <- read.table(file="~/Dropbox/Research/Papers/Ongoing/subtyping framework/Breast cancer classification/TCGA_RNA_seq/TCGA.BRCA.sampleMap_HiSeqV2", header = TRUE, row.names = 1)

## Subset both sets based on genes in common
testdata <- testdata[intersect(rownames(testdata), rownames(tcga_scaffold)),sample_names]
tcga_scaffold <- tcga_scaffold[intersect(rownames(testdata), rownames(tcga_scaffold)),]


# Simple spearman correlation to centroid (based on mean)
centroids <- do.call(cbind, lapply(unique(colnames(tcga_scaffold)), function(c) {rowMeans(tcga_scaffold[,colnames(tcga_scaffold)==c])}))
colnames(centroids) <- unique(colnames(tcga_scaffold))

cors <- do.call(rbind, lapply(1:10, function(i) {cor(testdata[,i], centroids, method = 'spearman')}))
rownames(cors) <- names(sample_names)
# heatmap(cors)

# cors <- do.call(rbind, lapply(1:10, function(i) {cor(testdata[,i], gtex_scaffold, method = 'spearman')}))
# rownames(cors) <- paste('sample', 1:10)

library(ComplexHeatmap)
# column_ha = HeatmapAnnotation(GTEx = class_sample)
tcga_plot <- Heatmap(cors, name = "Spearman's\ncorrelation", show_column_names = T, column_title_rot = 90, cluster_rows = F)
#              cluster_column_slices = T, top_annotation = column_ha, column_split = class_sample)

# pdf('../plot/Figure3_hm_tcga.pdf', height = 6, width = 8)
tcga_plot <- draw(tcga_plot, padding = unit(c(12.5, 2, 2, 2), "mm"), column_title = "TCGA heatmap")
# dev.off()
pdf('../figures_NOV23/Figure_3B.pdf', height = 8, width = 10)
draw(tcga_plot, padding = unit(c(43.7, 2, 2, 2), "mm"), column_title = "TCGA heatmap")
draw(tcga_plot, column_title = "GTEx heatmap")
dev.off()

# ggsave('../figures_NOV23/Figure_3B.png', height = 3.5, width = 6, dpi = 600)

# ggpubr::ggarrange(c(gtex_plot,  tcga_plot), nrow = 2,ncol=1, labels = c('A', 'B'))

################################################################################
################################### GTEX DGE ################################### 
################################################################################

#### DE GENES - THE SLACK LARS WAY ####
load(file="/home/people/s165827/data/Breast cancer classification/gtex.Rdata")

library(genefilter) # You may need to load this library for the coef function

# Create an empty dataframe to store log2fold changes
# log2fold_changes <- data.frame(row.names = rownames(tcga_scaffold))
log2fold_changes <- data.frame(row.names = rownames(gtex_scaffold))

# labels <- colnames(tcga_scaffold)
labels <- class_sample

# Find top DE probes (TPD)
groups <- list()
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
  significant_indices <- which(abs(log2fold_changes[, i]) > 0.2 & pvals < 0.05)
  # significant_indices <- which(pvals < 0.05)
  # top_10_significant_indices <- sort(pvals)[1:20]
  # groups[[i]] <- append(groups, as.vector(na.omit(names(pvals)[significant_indices])))
  groups[[i]] <- c(na.omit(names(pvals)[significant_indices]))
  # groups[[i]] <- names(top_10_significant_indices)

}

# TPD <- unique(features)
# tcga_scaffold_TPD <- tcga_scaffold[rownames(tcga_scaffold) %in% TPD, ]
# saveRDS(groups, 'rds/tcga_groups.rds')
# saveRDS(groups, 'rds/gtex_groups.rds')
# groups_tcga <- readRDS('rds/tcga_groups.rds')
groups_gtex <- readRDS('rds/gtex_groups.rds')
groups_gtex <- groups

# pca <- prcomp(t(tcga_scaffold_TPD), scale. = TRUE)

# De normaliserede ekspressionsværdier udregnes således:
gtex_scaffold_TPD_scaled <- NULL
for(i in 1:length(groups_gtex)) {
  pca <- prcomp(t(gtex_scaffold[groups_gtex[[i]],]), scale. = TRUE)
  gtex_scaffold_TPD_scaled <- rbind(gtex_scaffold_TPD_scaled, as.matrix(gtex_scaffold[groups_gtex[[i]],]/pca$sdev[1]^2))
}
gtex_scaffold_TPD_scaled

# saveRDS(tcga_scaffold_TPD_scaled, 'rds/tcga_scaffold_TPD_scaled.rds')
saveRDS(gtex_scaffold_TPD_scaled, 'rds/gtex_scaffold_TPD_scaled.rds')

# Rank??
gtex_scaffold_TPD_scaled_rank <- apply(gtex_scaffold_TPD_scaled, 2, rank)
# gtex_scaffold_TPD_scaled_rank <- gtex_scaffold_TPD_scaled

# PCA of tcga_scaffold_TPD_rank
pca <- prcomp(t(gtex_scaffold_TPD_scaled_rank), scale. = TRUE)
# pca_new <- scale(t(Bordet_RNA_core), pca$center, pca$scale) %*% pca$rotation
# df <- data.frame(PC1 = c(pca$x[,1], pca_new[,1]), PC2 = c(pca$x[,2], pca_new[,2]), source=c(rep("CITBCMST reference", ncol(CIT_full)), rep("Example data", ncol(Bordet_RNA_core))), sample=c(rep(NA,  ncol(CIT_full)), colnames(Bordet_RNA_core)))
df <- data.frame(PC1 = c(pca$x[,1]), PC2 = c(pca$x[,2]), tissue = class_sample)
var_exp <- summary(pca)$importance[2,1:2] * 100

# Plot
ggplot(df, aes(x = PC1, 
               y = PC2, 
               color = tissue)) +
  geom_point(size = 0.75) +
  ggtitle("TCGA PCA") +
  # ggrepel::geom_label_repel(size=1.75, show.legend = FALSE, alpha = 0.8, fontface = 'bold') +
  xlab(paste0('PC1 (', format(round(var_exp[1], 1), nsmall = 1), '% variance explained)')) + 
  ylab(paste0('PC2 (', format(round(var_exp[2], 1), nsmall = 1), '% variance explained)')) +
  theme_bw() 
  # theme(legend.position = "none")

################################################################################
################################ TCGA centroids ################################ 
################################################################################

# # Simple spearman correlation to centroid (based on mean)
centroids <- do.call(cbind, lapply(unique(colnames(tcga_scaffold_TPD_scaled_rank)), function(c) {rowMeans(tcga_scaffold_TPD_scaled_rank[,colnames(tcga_scaffold_TPD_scaled_rank)==c])}))
colnames(centroids) <- unique(colnames(tcga_scaffold_TPD_scaled))

# nc_pred <- c()
# for(i in 1:ncol(tcga_scaffold_TPD_scaled)) {
#   centroid_distances <- c()
#   for(ii in 1:length(centroids)) {
#     # centroid_distances <- c(centroid_distances, dist(rbind(CIT_TPD_rank[,i], centroids[[ii]])))
#     centroid_distances <- c(centroid_distances, DistancePair(tcga_scaffold_TPD_scaled[,i], centroids[[ii]]))
#   }
#   names(centroid_distances) <- names(centroids)
#   nc_pred <- c(nc_pred, names(which.min(centroid_distances)))
# }





################################################################################
################################### TCGA PCA ################################### 
################################################################################

# PCA of all
pca <- prcomp(t(tcga_scaffold), scale. = TRUE)
# pca_new <- scale(t(Bordet_RNA_core), pca$center, pca$scale) %*% pca$rotation
# df <- data.frame(PC1 = c(pca$x[,1], pca_new[,1]), PC2 = c(pca$x[,2], pca_new[,2]), source=c(rep("CITBCMST reference", ncol(CIT_full)), rep("Example data", ncol(Bordet_RNA_core))), sample=c(rep(NA,  ncol(CIT_full)), colnames(Bordet_RNA_core)))
df <- data.frame(PC1 = c(pca$x[,1]), PC2 = c(pca$x[,2]), tissue = colnames(tcga_scaffold))
var_exp <- summary(pca)$importance[2,1:2] * 100

# Plot
ggplot(df, aes(x = PC1, 
               y = PC2, 
               color = tissue)) +
  geom_point(size = 0.75) +
  ggtitle("TCGA PCA") +
  # ggrepel::geom_label_repel(size=1.75, show.legend = FALSE, alpha = 0.8, fontface = 'bold') +
  xlab(paste0('PC1 (', format(round(var_exp[1], 1), nsmall = 1), '% variance explained)')) + 
  ylab(paste0('PC2 (', format(round(var_exp[2], 1), nsmall = 1), '% variance explained)')) +
  theme_bw() +
  theme(legend.position = "none")








