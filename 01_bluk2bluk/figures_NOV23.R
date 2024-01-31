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
source('/home/people/s165827/data/Breast cancer classification/estimate_cbp.R')
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
################################### FIGURE 2 ################################### 
################################################################################
Bordet_RNA_tpm <- Bordet_RNA_tpm[,sample_names]

# Translate CIT and Bordet to gene names
library(genefilter)
library(hgu133plus2.db)
idx <- findLargest(gN = rownames(CIT_full), testStat = rowMedians(CIT_full), data = "hgu133plus2")
CIT_full_fl <- CIT_full[rownames(CIT_full) %in% idx, ]
rownames(CIT_full_fl) <- Bordet_annot[match(rownames(CIT_full_fl), Bordet_annot$Probe.Set.ID),]$Gene.Symbol
Bordet_RNA_fl <- Bordet_RNA_tpm[rownames(Bordet_RNA_tpm) %in% idx, ]
rownames(Bordet_RNA_fl) <- Bordet_annot[match(rownames(Bordet_RNA_fl), Bordet_annot$Probe.Set.ID),]$Gene.Symbol

# Run ESTIMATE on both
CIT_estimate <- ESTIMATE_fun(CIT_full_fl)
Bordet_estimate <- ESTIMATE_fun(Bordet_RNA_fl)
# range01 <- function(x){(x-min(CIT_estimate$TumorPurity))/(max(CIT_estimate$TumorPurity)-min(CIT_estimate$TumorPurity))}
# Bordet_estimate$TumorPurity <- range01(Bordet_estimate$TumorPurity)
# range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# CIT_estimate$TumorPurity <- range01(CIT_estimate$TumorPurity)

# Plot
cols = RColorBrewer::brewer.pal(10, 'Paired')[c(6:10,1:5)]
estimate_plot <- ggplot(CIT_estimate, aes(y=TumorPurity, x = 1)) +
  # geom_density() +
  geom_violin() + geom_boxplot(width = 0.1, outlier.shape = NA, coef = 0) + coord_flip(ylim = c(0.03,0.97)) +
  ggtitle("Purity distribution for CITBCMST set + example data") +
  geom_point(data = Bordet_estimate, aes(y = TumorPurity, x = 1 + rnorm(10, sd = 0.05))) +
  ggrepel::geom_text_repel(data = Bordet_estimate, aes(y = TumorPurity), label = paste0('sample ', 1:10), max.overlaps = 10, min.segment.length = 0.5, force = 30, direction = 'y') +
  # geom_vline(xintercept=Bordet_estimate$TumorPurity, linetype="dashed") +
  # annotate("text", label = paste0("sample ", 1:10), x=seq(0.1,1,0.1), y=Bordet_estimate$TumorPurity) +
  ylab('Tumor purity') + xlab('') + scale_x_discrete(labels = NULL, breaks = NULL) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) 
  
estimate_plot

ggsave(estimate_plot, filename = "../figures_NOV23/Figure_2.png", height = 3.5, width = 6, dpi = 600)



################################################################################
################################### FIGURE 3 ################################### 
################################################################################

## A) GTEx

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

## B) TCGA

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
################################### FIGURE 4 ################################### 
################################################################################

# Compare Bordet RNA seq with CIT
# load(file="Bordet_paired_data/Bordet.rdata")
# load("CITBCMST.rdata")

Bordet_RNA_core <- Bordet_RNA_tpm[rownames(Bordet_RNA_tpm) %in% rownames(CIT_core),sample_names]
Bordet_RNA_core <- Bordet_RNA_core[match(rownames(CIT_core), rownames(Bordet_RNA_core)),]
colnames(Bordet_RNA_core) <- paste("sample", 1:10)

## A) not corrected
# PCA of all
pca <- prcomp(t(CIT_core), scale. = TRUE)
pca_new <- scale(t(Bordet_RNA_core), pca$center, pca$scale) %*% pca$rotation
df <- data.frame(PC1 = c(pca$x[,1], pca_new[,1]), PC2 = c(pca$x[,2], pca_new[,2]), source=c(rep("CITBCMST reference", ncol(CIT_full)), rep("Example data", ncol(Bordet_RNA_core))), sample=c(rep(NA,  ncol(CIT_full)), colnames(Bordet_RNA_core)))
var_exp <- summary(pca)$importance[2,1:2] * 100

# Plot
library(ggplot2)
uncorrected <- ggplot(df, aes(x=PC1, y=PC2, color=source, label=sample)) +
  geom_point(size = 0.75) +
  ggtitle("Uncorrected") +
  scale_color_manual(values = c('#f8766d', '#1f7c81')) + 
  ggrepel::geom_label_repel(size=1.75, show.legend = FALSE, alpha = 0.8, fontface = 'bold') +
  xlab(paste0('PC1 (', format(round(var_exp[1], 1), nsmall = 1), '% variance explained)')) + ylab(paste0('PC2 (', format(round(var_exp[2], 1), nsmall = 1), '% variance explained)')) +
  theme_bw() +
  theme(legend.position = "none") 


## B) COMBAT corrected
library(sva)
library(preprocessCore)
# Bordet_RNA_qnorm <- normalize.quantiles.use.target(Bordet_RNA_tpm[,sample_names], target = rowMeans(CIT_full))
# rownames(Bordet_RNA_qnorm) <- rownames(Bordet_RNA_tpm); colnames(Bordet_RNA_qnorm) <- paste("sample", 1:10)
# exprs_bc <- ComBat(dat = cbind(CIT_full, Bordet_RNA_qnorm), batch = c(rep(1, ncol(CIT_full)), rep(2, ncol(Bordet_RNA_qnorm))), ref.batch = 1)
# CIT_core_bc <- exprs_bc[rownames(exprs_bc) %in% rownames(CIT_core),colnames(exprs_bc) %in% colnames(CIT_core)]
# Bordet_RNA_core_bc <- exprs_bc[rownames(exprs_bc) %in% rownames(CIT_core),colnames(exprs_bc) %in% colnames(Bordet_RNA_qnorm)]

# PCA of all
pca <- prcomp(t(CIT_core_bc), scale. = TRUE)
pca_new <- scale(t(Bordet_RNA_core_bc), pca$center, pca$scale) %*% pca$rotation
df <- data.frame(PC1 = c(pca$x[,1], pca_new[,1]), PC2 = c(pca$x[,2], pca_new[,2]), Source=c(rep("CITBCMST reference", ncol(CIT_core_bc)), rep("Example data", ncol(Bordet_RNA_core_bc))), sample=c(rep(NA,  ncol(CIT_core_bc)), colnames(Bordet_RNA_core_bc)))
var_exp <- summary(pca)$importance[2,1:2] * 100

# Plot
combat <- ggplot(df, aes(x=PC1, y=PC2, color=Source, label=sample)) +
  geom_point(size = 0.75) +
  ggtitle("ComBat") +
  scale_color_manual(values = c('#f8766d', '#1f7c81')) + 
  ggrepel::geom_label_repel(size=1.75, show.legend = FALSE, alpha = 0.8, fontface = 'bold') +
  xlab(paste0('PC1 (', format(round(var_exp[1], 1), nsmall = 1), '% variance explained)')) + ylab(paste0('PC2 (', format(round(var_exp[2], 1), nsmall = 1), '% variance explained)')) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "vertical")


## C) Rank transformed
Bordet_RNA_rank <- apply(Bordet_RNA_core, 2, rank)
CIT_rank <- apply(CIT_core, 2, rank)

# PCA of both
pca <- prcomp(t(CIT_rank), scale. = TRUE)
pca_new <- scale(t(Bordet_RNA_rank), pca$center, pca$scale) %*% pca$rotation
df <- data.frame(PC1 = c(pca$x[,1], pca_new[,1]), PC2 = c(pca$x[,2], pca_new[,2]), Source=c(rep("CITBCMST reference", ncol(CIT_rank)), rep("Example data", ncol(Bordet_RNA_rank))), sample=c(rep(NA,  ncol(CIT_rank)), colnames(Bordet_RNA_rank)))
var_exp <- summary(pca)$importance[2,1:2] * 100

# Plot
ranktransformed <- ggplot(df, aes(x=PC1, y=PC2, color=Source, label=sample)) +
  geom_point(size = 0.75) +
  ggtitle("Rank") +
  scale_color_manual(values = c('#f8766d', '#1f7c81')) + 
  xlab(paste0('PC1 (', format(round(var_exp[1], 1), nsmall = 1), '% variance explained)')) + ylab(paste0('PC2 (', format(round(var_exp[2], 1), nsmall = 1), '% variance explained)')) +
  ggrepel::geom_label_repel(size=1.75, show.legend = FALSE, alpha = 0.8, fontface = 'bold') +
  theme_bw() + 
  theme(legend.position = 'none')


## Plot all
# cowplot::plot_grid(uncorrected, combat, ranktransformed, labels = c("A", "B", "C"), ncol=3)
plot_be <- uncorrected + combat + ranktransformed +
  plot_layout(design = 'ABC') +
  plot_annotation(tag_levels = list(c('A', 'B', 'C'))) &
  theme(plot.tag = element_text(size = 8, face = "bold"), #16
        plot.title = element_text(size = 8),
        axis.title.x = element_text(size = 7),
        axis.text.x = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 5),
        legend.title = element_text(size=8), 
        legend.text = element_text(size=7)
        )  
plot_be
# TOO BIG 
# ggsave(plot_be, filename = '../figures_NOV23/Figure_4.png', height = 5, width = 12, dpi = 600)
ggsave(plot_be, filename = '../figures_NOV23/Figure_4.png', height = 3.5, width = 7, dpi = 600)


################################################################################
################################### FIGURE 5 ################################### 
################################################################################

## Test CIT subtype copy vs whole space Combat vs rank vs biological process transformation
library(sva)

# load("~/Dropbox/Breast cancer classification/CITBCMST.rdata")
lumb_sub <- CIT_full[,CIT_classes == "lumB"]

exprs_bc <- ComBat(dat = cbind(CIT_full, lumb_sub), batch = c(rep(1, ncol(CIT_full)), rep(2, ncol(lumb_sub))), ref.batch = 1)
CIT_core_bc <- exprs_bc[rownames(exprs_bc) %in% rownames(CIT_core),1:ncol(CIT_full)]
lumb_core_bc <- exprs_bc[rownames(exprs_bc) %in% rownames(CIT_core),(ncol(CIT_full)+1):ncol(exprs_bc)]

pca <- prcomp(t(CIT_core_bc), scale. = TRUE)
pca_new <- scale(t(lumb_core_bc), pca$center, pca$scale) %*% pca$rotation
var_exp <- summary(pca)$importance[2,1:2] * 100
df <- data.frame(PC1 = c(pca$x[,1], pca_new[,1]), PC2 = c(pca$x[,2], pca_new[,2]), Source=c(rep("Full CITBCMST", ncol(CIT_full)), rep("lumB only", ncol(lumb_core_bc))))
df$Subtype <- factor(c(CIT_classes, rep("lumB", ncol(lumb_core_bc))), levels = cit_lvls)
df$sample <- c(colnames(CIT_core_bc), colnames(lumb_core_bc))

# Plot
combat <- ggplot(df, aes(x=PC1, y=PC2, color=Subtype, shape=Source)) +
  geom_point(size = 1) +
  scale_shape_manual(values = c(16,4)) + scale_color_manual(values = cit_colors) +
  xlab(paste0('PC1 (', format(round(var_exp[1], 1), nsmall = 1), '% variance explained)')) + ylab(paste0('PC2 (', format(round(var_exp[2], 1), nsmall = 1), '% variance explained)')) +
  geom_line(aes(group = sample), alpha = 0.4) +
  ggtitle("ComBat") +
  theme_bw()

## Rank transformed
CIT_rank <- apply(CIT_full[rownames(CIT_full) %in% rownames(CIT_core),], 2, rank)
CIT_lumc_rank <- apply(lumb_sub[rownames(lumb_sub) %in% rownames(CIT_core),], 2, rank)

# PCA of both
pca <- prcomp(t(CIT_rank), scale. = TRUE)
pca_new <- scale(t(CIT_lumc_rank), pca$center, pca$scale) %*% pca$rotation
var_exp <- summary(pca)$importance[2,1:2] * 100
df <- data.frame(PC1 = c(pca$x[,1], pca_new[,1]), PC2 = c(pca$x[,2], pca_new[,2]), Source=c(rep("Full CITBCMST", ncol(CIT_rank)), rep("lumB only", ncol(CIT_lumc_rank))))
df$Subtype <- factor(c(CIT_classes, rep("lumB", ncol(CIT_lumc_rank))), levels = cit_lvls)
df$sample <- c(colnames(CIT_rank), colnames(CIT_lumc_rank))

# Plot
rank <- ggplot(df, aes(x=PC1, y=PC2, color=Subtype, shape=Source)) +
  geom_point(size = 1) +
  scale_shape_manual(values = c(16,4)) + scale_color_manual(values = cit_colors) +
  xlab(paste0('PC1 (', format(round(var_exp[1], 1), nsmall = 1), '% variance explained)')) + ylab(paste0('PC2 (', format(round(var_exp[2], 1), nsmall = 1), '% variance explained)')) +
  geom_line(aes(group = sample), alpha = 0.4) +
  ggtitle("Rank") +
  theme_bw()

# cowplot::plot_grid(combat, rank, labels = c("A", "B"))
plot_bc <- combat + rank +
  plot_layout(design = 'AB', guides = 'collect') +
  plot_annotation(tag_levels = list(c('A', 'B'))) &
  theme(plot.tag = element_text(size = 10, face = "bold"), #16
        plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        axis.text.y = element_text(size = 6),
        legend.title=element_text(size=8), 
        legend.text =element_text(size=7)
        )

# ggsave(plot_bc, filename = '../figures_NOV23/Figure_5.png', height = 5, width = 11, dpi = 600)
ggsave(plot_bc, filename = '../figures_NOV23/Figure_5.png', height = 3.2, width = 7, dpi = 600)

################################################################################
################################### FIGURE 6 ################################### 
################################################################################

source("performance_eval.R")
# 
# # ## loo distance to centroid on TDP
# centroids <- list()
# for (i in unique(CIT_classes)) {
#   centroids[[i]] <- rowMeans(CIT_TPD[,CIT_classes==i])
#   
# }
# centroids_rank <- lapply(centroids, rank)
# 
# library(rankdist)
# nc_pred <- c()
# for(i in 1:ncol(CIT_TPD_rank)) {
#   centroid_distances <- c()
#   for(ii in 1:length(centroids)) {
#     # centroid_distances <- c(centroid_distances, dist(rbind(CIT_TPD_rank[,i], centroids[[ii]])))
#     centroid_distances <- c(centroid_distances, DistancePair(CIT_TPD_rank[,i], centroids[[ii]]))
#   }
#   names(centroid_distances) <- names(centroids)
#   nc_pred <- c(nc_pred, names(which.min(centroid_distances)))
# }
# 
# perf <- performance_eval(CIT_classes, nc_pred)
# perf[[1]]
# perf[[2]]

### FIGURE 6: fuzzy classification of Bordet samples using best performing method from above
load("~/Dropbox/Breast cancer classification/Bordet_paired_data/Bordet.Rdata")
load("~/Dropbox/Breast cancer classification/CITBCMST.rdata")

## Predict
features <- c()
for(i in unique(CIT_classes)) {
  pvals <- apply(CIT_core, 1, function(x) wilcox.test(x[CIT_classes == i], x[!CIT_classes == i])$p.value)
  pvals <- p.adjust(pvals)
  features <- append(features, as.vector(na.omit(names(pvals)[pvals<0.05])))
  print(i)
}
TPD <- unique(features)
CIT_TPD <- CIT_full[rownames(CIT_full) %in% TPD, ]

## Rank
Bordet_RNA_core <- Bordet_RNA_tpm[rownames(Bordet_RNA_tpm) %in% TPD,sample_names]
Bordet_RNA_rank <- apply(Bordet_RNA_core, 2, rank)
colnames(Bordet_RNA_rank) <- paste0("sample ", 1:10)
CIT_TPD_rank <- apply(CIT_TPD, 2, rank)

centroids <- list()
for (i in unique(CIT_classes)) {
  centroids[[i]] <- rowMeans(CIT_TPD[,CIT_classes==i])
  
}
centroids_rank <- lapply(centroids, rank)

library(rankdist)
cd <- c()
nc_pred <- c()
for(i in 1:ncol(Bordet_RNA_rank)) {
  centroid_distances <- c()
  for(ii in 1:length(centroids_rank)) {
    # centroid_distances <- c(centroid_distances, dist(rbind(CIT_TPD_rank[,i], centroids[[ii]])))
    centroid_distances <- c(centroid_distances, DistancePair(Bordet_RNA_rank[,i], centroids_rank[[ii]]))
  }
  # names(centroid_distances) <- names(centroids)
  cd <- rbind(cd, centroid_distances)
  nc_pred <- c(nc_pred,names(centroids_rank)[which.min(centroid_distances)])
}
rownames(cd) <- colnames(Bordet_RNA_rank)
colnames(cd) <- names(centroids_rank)


## Plot PCA
pca <- prcomp(t(CIT_TPD_rank), scale. = TRUE)
var_exp <- summary(pca)$importance[2,1:2] * 100
pca_new <- scale(t(Bordet_RNA_rank), pca$center, pca$scale) %*% pca$rotation
df <- data.frame(PC1 = c(pca$x[,1], pca_new[,1]), PC2 = c(pca$x[,2], pca_new[,2]), Source=c(rep("CITBCMST", ncol(CIT_TPD_rank)), rep("Example data", ncol(Bordet_RNA_rank))), Subtype=factor(c(CIT_classes, nc_pred), levels = cit_lvls))
df$sample <- rownames(df); df[df$Source=="CITBCMST",]$sample <- NA

# Plot
pca <- ggplot(df, aes(x=PC1, y=PC2, color=Subtype, shape=Source, label=sample, size=Source)) +
  geom_point(size = 1) +
  scale_shape_manual(values = c(16,4)) + scale_size_manual(values = c(2,3.5)) + scale_color_manual(values = cit_colors) +
  xlab(paste0('PC1 (', format(round(var_exp[1], 1), nsmall = 1), '% variance explained)')) + ylab(paste0('PC2 (', format(round(var_exp[2], 1), nsmall = 1), '% variance explained)')) +
  ggrepel::geom_label_repel(size=2, show.legend = FALSE, alpha = 0.8, fontface = 'bold') +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.box="vertical",
        legend.margin=margin())
pca

## Plot heatmap
library(reshape2)
cds <- t(apply(cd, 1, scale))
colnames(cds) <- colnames(cd)
df <- melt(cds)
colnames(df) <- c("Sample", "Subtype", "distance")
df$Subtype <- factor(df$Subtype, levels = cit_lvls)
df$Sample <- factor(str_split_i(df$Sample, ' ', 2), levels = c('10', '9', '8', '7', '6', '5', '4', '3', '2', '1'))

heat <- ggplot(df, aes(Subtype, Sample, fill= distance)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient2(low="dark red", mid="white", high="dark blue", midpoint = mean(df$distance)) +
  labs(fill = 'Scaled\ndistance') +
  geom_rect(size=0.3, fill=NA, colour="black", xmin = 5.5, xmax = 6.5, ymin = 0.5, ymax = 1.5) + #10
  geom_rect(size=0.3, fill=NA, colour="black", xmin = 4.5, xmax = 5.5, ymin = 1.5, ymax = 2.5) + # 9
  geom_rect(size=0.3, fill=NA, colour="black", xmin = 3.5, xmax = 4.5, ymin = 2.5, ymax = 3.5) + # 8
  geom_rect(size=0.3, fill=NA, colour="black", xmin = 1.5, xmax = 2.5, ymin = 3.5, ymax = 4.5) + # 7
  geom_rect(size=0.3, fill=NA, colour="black", xmin = 0.5, xmax = 1.5, ymin = 4.5, ymax = 5.5) + # 6
  geom_rect(size=0.3, fill=NA, colour="black", xmin = 0.5, xmax = 1.5, ymin = 5.5, ymax = 6.5) + # 5
  geom_rect(size=0.3, fill=NA, colour="black", xmin = 1.5, xmax = 2.5, ymin = 6.5, ymax = 7.5) + # 4
  geom_rect(size=0.3, fill=NA, colour="black", xmin = 0.5, xmax = 1.5, ymin = 7.5, ymax = 8.5) + # 3
  geom_rect(size=0.3, fill=NA, colour="black", xmin = 4.5, xmax = 5.5, ymin = 8.5, ymax = 9.5) + # 2
  geom_rect(size=0.3, fill=NA, colour="black", xmin = 2.5, xmax = 3.5, ymin = 9.5, ymax = 10.5) + # 1
  theme_bw(base_size=8) +
  theme(
    #bold font for legend text
    legend.text=element_text(face="bold"),
    #set thickness of axis ticks
    axis.ticks=element_line(size=0.4),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank(),
    axis.title.x = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 6),
    legend.title=element_text(size=7), 
    legend.position = 'bottom'
    )

# cowplot::plot_grid(pca, heat, labels = c("A", "B"))
plot_fuzzy <- pca + heat +
  plot_layout(design = 'AB') +
  plot_annotation(tag_levels = list(c('A', 'B'))) &
  theme(plot.tag = element_text(size = 10, face = "bold"), #16
        plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        axis.text.y = element_text(size = 6),
        legend.title=element_text(size=7), 
        legend.text =element_text(size=6)
  )
  # theme(
  #   plot.tag = element_text(size = 8, face = "bold"), #16
  #   plot.title = element_text(size = 8),
  #   axis.title.x = element_text(size = 8),
  #   axis.text.x = element_text(size = 6),
  #   axis.title.y = element_text(size = 8),
  #   axis.text.y = element_text(size = 6),
  #   legend.title=element_text(size=8), 
  #   legend.text =element_text(size=6)
  #       )
plot_fuzzy
# ggsave(plot_fuzzy, filename = '../figures_NOV23/Figure_6.png', height = 5, width = 12, dpi = 600)
ggsave(plot_fuzzy, filename = '../figures_NOV23/Figure_6.png', height = 4, width = 6, dpi = 600)

################################################################################
################################### FIGURE 7 ################################### 
################################################################################

## A) PCA of TCGA space w Bordet samples projected
# load("brca.Rdata")
# load("PAM50genes.Rdata")
brca <- brca[rownames(brca) %in% pam50_genes, colnames(brca) %in% brca_pam50$sample]
brca_pam50 <- brca_pam50[brca_pam50$sample %in% colnames(brca),]
subtypes <- brca_pam50[match(colnames(brca), brca_pam50$sample),]$PAM50

TdeG <- pam50_genes

# Sub and rank
brca_rank <- apply(brca[rownames(brca) %in% TdeG,], 2, rank)

# load("Bordet_RNA_gene_tpm.Rdata")
Bordet_array_TdeG <- Bordet_array[rownames(Bordet_array) %in% TdeG, sample_names]
Bordet_RNA_TdeG <- Bordet_RNA_gene_tpm[rownames(Bordet_RNA_gene_tpm) %in% TdeG, sample_names] # RNA
colnames(Bordet_RNA_TdeG) <- paste('sample', 1:10)

Bordet_RNA_TdeG_rank <- apply(Bordet_RNA_TdeG, 2, rank)
brca_TdeG_rank <- apply(brca[rownames(brca) %in% rownames(Bordet_RNA_TdeG_rank), ], 2, rank)
Bordet_RNA_TdeG_rank <- Bordet_RNA_TdeG_rank[match(rownames(brca_TdeG_rank), rownames(Bordet_RNA_TdeG_rank)),]

Bordet_array_TdeG_rank <- apply(Bordet_RNA_TdeG, 2, rank)
# brca_TdeG_rank <- apply(brca[rownames(brca) %in% rownames(Bordet_array_TdeG), ], 2, rank)
Bordet_array_TdeG_rank <- Bordet_array_TdeG_rank[match(rownames(brca_TdeG_rank), rownames(Bordet_array_TdeG_rank)),]

# Classify Bordet
centroids <- list()
for (i in unique(subtypes)) {
  centroids[[i]] <- rowMeans(brca[rownames(brca) %in% rownames(Bordet_RNA_TdeG_rank), subtypes==i])
  centroids[[i]] <- apply(brca[rownames(brca) %in% TdeG, subtypes==i], 1, median)
  
}
centroids_rank <- lapply(centroids, rank)

#### ND PRED RNA
library(rankdist)
nc_pred_RNA <- c()
cd <- NULL
for(i in 1:ncol(Bordet_RNA_TdeG_rank)) {
  centroid_distances <- c()
  for(ii in 1:length(centroids_rank)) {
    # centroid_distances <- c(centroid_distances, dist(rbind(brca_TdeG_rank[,i], centroids_rank[[ii]]))) # nc_pred_B
    centroid_distances <- c(centroid_distances, DistancePair(Bordet_RNA_TdeG_rank[,i], centroids_rank[[ii]])) #nc_pred_A
  }
  names(centroid_distances) <- names(centroids)
  nc_pred_RNA <- c(nc_pred_RNA, names(which.min(centroid_distances)))
  cd <- rbind(cd, centroid_distances)
}
rownames(cd) <- colnames(Bordet_RNA_TdeG_rank)
colnames(cd) <- names(centroids_rank)

#### ND PRED ARAAY
nc_pred_array <- c()
cd <- NULL
for(i in 1:ncol(Bordet_array_TdeG_rank)) {
  centroid_distances <- c()
  for(ii in 1:length(centroids_rank)) {
    # centroid_distances <- c(centroid_distances, dist(rbind(brca_TdeG_rank[,i], centroids_rank[[ii]]))) # nc_pred_B
    centroid_distances <- c(centroid_distances, DistancePair(Bordet_array_TdeG_rank[,i], centroids_rank[[ii]])) #nc_pred_A
  }
  names(centroid_distances) <- names(centroids)
  nc_pred_array <- c(nc_pred_array, names(which.min(centroid_distances)))
  cd <- rbind(cd, centroid_distances)
}
rownames(cd) <- colnames(Bordet_array_TdeG_rank)
colnames(cd) <- names(centroids_rank)

table(nc_pred_array == nc_pred_RNA)



# Plot PCA
## Plot PCA
pca <- prcomp(t(brca_TdeG_rank), scale. = TRUE)
pca_new <- scale(t(Bordet_RNA_TdeG_rank), pca$center, pca$scale) %*% pca$rotation
df <- data.frame(PC1 = c(pca$x[,1], pca_new[,1]), PC2 = c(pca$x[,2], pca_new[,2]), Source=c(rep("TCGA (reference)", ncol(brca_TdeG_rank)), rep("Example data", ncol(Bordet_RNA_TdeG_rank))), Subtype=factor(c(subtypes, nc_pred_RNA), levels = pam50_lvls)) #nc_pred_RNA == nc_pred_array for the 10 samples
df$sample <- rownames(df); df[df$Source=="TCGA (reference)",]$sample <- NA
df$Source <- factor(df$Source, levels = c('TCGA (reference)', 'Example data'))
var_exp <- summary(pca)$importance[2,1:2] * 100

# Plot
pca <- ggplot(df, aes(x=PC1, y=PC2, color=Subtype, shape=Source, label=sample, size=Source)) +
  geom_point(size = 1) +
  scale_shape_manual(values = c(16,4)) + scale_size_manual(values = c(2,3.5)) + scale_color_manual(values = pam50_colors) +
  xlab(paste0('PC1 (', format(round(var_exp[1], 1), nsmall = 1), '% variance explained)')) + ylab(paste0('PC2 (', format(round(var_exp[2], 1), nsmall = 1), '% variance explained)')) +
  ggrepel::geom_label_repel(size=2, show.legend = FALSE, alpha = 0.8, fontface = 'bold') +
  guides(color = guide_legend(order = 1, nrow = 2),
         shape = guide_legend(order = 2),
         size = guide_legend(order = 2)) +
  theme_bw() + 
  theme(legend.position = 'bottom',
        legend.box="vertical",
        legend.margin=margin()) 

pca
# ggsave('../helenewegener/Desktop/subtyping/analysis/plot/PAM50_pca.pdf')

# Plot heat
library(reshape2)
cds <- t(apply(cd, 1, scale))
colnames(cds) <- colnames(cd)
df <- melt(cds)
colnames(df) <- c("Sample", "Subtype", "distance")
df$Subtype <- factor(df$Subtype, levels = pam50_lvls)
df$Sample <- factor(str_split_i(df$Sample, ' ', 2), levels = c('10', '9', '8', '7', '6', '5', '4', '3', '2', '1'))

heat <- ggplot(df, aes(Subtype, Sample, fill= distance)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient2(low="dark red", mid="white", high="dark blue", midpoint = mean(df$distance)) +
  labs(fill = 'Scaled\ndistance') +
  geom_rect(size=0.3, fill=NA, colour="black", xmin = 4.5, xmax = 5.5, ymin = 0.5, ymax = 1.5) + #10
  geom_rect(size=0.3, fill=NA, colour="black", xmin = 0.5, xmax = 1.5, ymin = 1.5, ymax = 2.5) + # 9
  geom_rect(size=0.3, fill=NA, colour="black", xmin = 2.5, xmax = 3.5, ymin = 2.5, ymax = 3.5) + # 8
  geom_rect(size=0.3, fill=NA, colour="black", xmin = 1.5, xmax = 2.5, ymin = 3.5, ymax = 4.5) + # 7
  geom_rect(size=0.3, fill=NA, colour="black", xmin = 1.5, xmax = 2.5, ymin = 4.5, ymax = 5.5) + # 6
  geom_rect(size=0.3, fill=NA, colour="black", xmin = 1.5, xmax = 2.5, ymin = 5.5, ymax = 6.5) + # 5
  geom_rect(size=0.3, fill=NA, colour="black", xmin = 1.5, xmax = 2.5, ymin = 6.5, ymax = 7.5) + # 4
  geom_rect(size=0.3, fill=NA, colour="black", xmin = 1.5, xmax = 2.5, ymin = 7.5, ymax = 8.5) + # 3
  geom_rect(size=0.3, fill=NA, colour="black", xmin = 3.5, xmax = 4.5, ymin = 8.5, ymax = 9.5) + # 2
  geom_rect(size=0.3, fill=NA, colour="black", xmin = 1.5, xmax = 2.5, ymin = 9.5, ymax = 10.5) + # 1
  theme_bw(base_size=8) +
  theme(
    #bold font for legend text
    legend.text=element_text(face="bold"),
    #set thickness of axis ticks
    axis.ticks=element_line(size=0.4),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank(),
    axis.title.x = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 6),
    legend.title=element_text(size=7),
    legend.position = 'bottom'#,
    # legend.box="vertical"
    )

# cowplot::plot_grid(pca, heat, labels = c("A", "B"))
plot_fuzzy_pam50 <- pca + heat +
  plot_layout(design = 'AB') +
  plot_annotation(tag_levels = list(c('A', 'B'))) &
  theme(
    plot.tag = element_text(size = 8, face = "bold"), #16
    plot.title = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 6),
    legend.title=element_text(size=7), 
    legend.text =element_text(size=6)
  )
plot_fuzzy_pam50

# ggsave(plot_fuzzy_pam50, filename = '../figures_NOV23/Figure_7.png', height = 5, width = 12, dpi = 600)
ggsave(plot_fuzzy_pam50, filename = '../figures_NOV23/Figure_7.png', height = 4, width = 6, dpi = 600)

################################################################################
################################### FIGURE 8 ################################### 
################################################################################

load("Bordet.rdata")
load("CITBCMST.rdata")

library(genefilter)
library(hgu133plus2.db)
idx <- findLargest(gN = rownames(CIT_full), testStat = rowMedians(CIT_full), data = "hgu133plus2")
CIT_full_fl <- CIT_full[rownames(CIT_full) %in% idx, ]
rownames(CIT_full_fl) <- Bordet_annot[match(rownames(CIT_full_fl), Bordet_annot$Probe.Set.ID),]$Gene.Symbol
CIT_full_rank <- apply(CIT_full_fl, 2, rank)
Bordet_RNA_fl <- Bordet_RNA_tpm[rownames(Bordet_RNA_tpm) %in% idx, sample_names]
rownames(Bordet_RNA_fl) <- Bordet_annot[match(rownames(Bordet_RNA_fl), Bordet_annot$Probe.Set.ID),]$Gene.Symbol
# colnames(Bordet_RNA_fl) <- paste0("sample ", 1:10)
Bordet_RNA_rank <- apply(Bordet_RNA_fl, 2, rank)

library(ggExtra)


# Density plots
df <- data.frame(HER2 = CIT_full_rank[rownames(CIT_full_rank) == "ERBB2",], sample = colnames(CIT_full_rank), source=rep("CIT", ncol(CIT_full_rank)))
HER2 <- ggplot(df, aes(x=HER2)) +
  geom_density() +
  xlab(bquote(italic(.('ERBB2'))~'expression rank')) +
  geom_vline(xintercept=Bordet_RNA_rank[rownames(Bordet_RNA_rank) == "ERBB2",], linetype="dashed") +
  annotate("text", label = paste('sample', 1:10), y=seq(0, 3.9e-04, by = 3.9e-04/10)[-1], x=Bordet_RNA_rank[rownames(Bordet_RNA_rank) == "ERBB2",],
           size = 3.5) +
  ggtitle(bquote('Expression rank of'~italic(.('ERBB2'))~'gene')) + ylab('Density') +
  theme_bw()


df <- data.frame(ER = CIT_full_rank[rownames(CIT_full_rank) == "ESR1",], sample = colnames(CIT_full_rank), source=rep("CIT", ncol(CIT_full_rank)))
ER <- ggplot(df, aes(x=ER)) +
  geom_density() +
  xlab(bquote(italic(.('ESR1'))~'expression rank')) +
  geom_vline(xintercept=Bordet_RNA_rank[rownames(Bordet_RNA_rank) == "ESR1",], linetype="dashed") +
  annotate("text", label = paste('sample', 1:10), y=seq(0, 3.7e-04, by = 3.7e-04/10)[-1], x=Bordet_RNA_rank[rownames(Bordet_RNA_rank) == "ESR1",], 
           size = 3.5) +
  ggtitle(bquote('Expression rank of'~italic(.('ESR1'))~'gene')) + ylab('Density') +
  theme_bw()

df <- data.frame(PR = CIT_full_rank[rownames(CIT_full_rank) == "PGR",], sample = colnames(CIT_full_rank), source=rep("CIT", ncol(CIT_full_rank)))
PR <- ggplot(df, aes(x=PR)) +
  geom_density() +
  xlab(bquote(italic(.('PGR'))~'expression rank')) +
  geom_vline(xintercept=Bordet_RNA_rank[rownames(Bordet_RNA_rank) == "PGR",], linetype="dashed") +
  annotate("text", label = paste('sample', 1:10), y=seq(0, 9e-05, by = 9e-05/10)[-1], x=Bordet_RNA_rank[rownames(Bordet_RNA_rank) == "PGR",],
           size = 3.5) +
  ggtitle(bquote('Expression rank of'~italic(.('PGR'))~'gene')) + ylab('Density') +
  theme_bw()

# cowplot::plot_grid(HER2, ER, PR, labels = c("A", "B", "C"), ncol=1)
plot_den_rec <- HER2 + ER + PR +
  plot_layout(design = 'A\nB\nC') +
  plot_annotation(tag_levels = list(c('A', 'B', 'C'))) &
  theme(plot.tag = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 12))
plot_den_rec

# ggsave(plot_den_rec, filename = '../figures_NOV23/Figure_8.png', height = 10, width = 9, dpi = 600)
ggsave(plot_den_rec, filename = '../figures_NOV23/Figure_8.png', height = 7.8, width = 7, dpi = 600)












