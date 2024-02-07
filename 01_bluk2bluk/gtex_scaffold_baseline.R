library(spaceRATScaffolds)
library(spaceRAT)
library(tidyverse)

load(file="/home/people/s165827/data/Breast cancer classification/gtex.Rdata")
gtex_pheno <- as.data.frame(class_sample)
rownames(gtex_pheno) <- colnames(gtex_scaffold)

################################################################################
################################## Baseline ################################### 
################################################################################

gtex_Scaffold_baseline <- buildScaffold(gtex_scaffold, 
                                        pheno = gtex_pheno,
                                        colname = 'class_sample',
                                        data = 'exprs', 
                                        subset_deg = TRUE,
                                        threshold = 10,
                                        pca_scale = TRUE)

plotScaffold(gtex_Scaffold_baseline, dims = c(1,2), title = 'gtex scaffold baseline') + theme(legend.position = "none")
ggsave('plots/gtex_scaffold_baseline_dim12.png', width = 10, height = 8.72)
plotScaffold(gtex_Scaffold_baseline, dims = c(2,3), title = 'gtex scaffold baseline') + theme(legend.position = "none")
ggsave('plots/gtex_scaffold_baseline_dim23.png', width = 10, height = 8.72)
plotScaffold(gtex_Scaffold_baseline, dims = c(1,3), title = 'gtex scaffold baseline') + theme(legend.position = "none")
ggsave('plots/gtex_scaffold_baseline_dim13.png', width = 10, height = 8.72)

################################################################################
########################### Remove Kidney Outlier ############################## 
################################################################################

# Identify outlier 
# plotScaffold(gtex_scaffold_baseline, dims = c(1,2), title = 'gtex scaffold baseline', plot_mode = 'tiny_label') + theme(legend.position = "none")
outlier <- gtex_Scaffold_baseline$pca$x[, 'PC1'] %>% sort() %>% head(1) %>% names()

# Remove outlier
gtex_scaffold_wo_outlier <- gtex_scaffold[, colnames(gtex_scaffold) != outlier]
gtex_pheno_wo_outlier <- as.data.frame(gtex_pheno[rownames(gtex_pheno) != outlier ,])
rownames(gtex_pheno_wo_outlier) <- rownames(gtex_pheno)[rownames(gtex_pheno) != outlier]
colnames(gtex_pheno_wo_outlier) <- c('class_sample')

saveRDS(gtex_scaffold_wo_outlier, 'rds/gtex_scaffold_wo_outlier.rds')
saveRDS(gtex_pheno_wo_outlier, 'rds/gtex_pheno_wo_outlier.rds')

gtex_Scaffold_wo_outlier <- buildScaffold(gtex_scaffold_wo_outlier, 
                                          pheno = gtex_pheno_wo_outlier,
                                          colname = 'class_sample',
                                          data = 'exprs', 
                                          subset_deg = TRUE,
                                          threshold = 10,
                                          pca_scale = TRUE)

plotScaffold(gtex_Scaffold_wo_outlier, dims = c(1,2), title = 'gtex scaffold w/o outlier') + theme(legend.position = "none")
ggsave('plots/gtex_scaffold_wo_outlier_dim12.png', width = 10, height = 8.72)
plotScaffold(gtex_Scaffold_wo_outlier, dims = c(2,3), title = 'gtex scaffold w/o outlier') + theme(legend.position = "none")
ggsave('plots/gtex_scaffold_wo_outlier_dim23.png', width = 10, height = 8.72)
plotScaffold(gtex_Scaffold_wo_outlier, dims = c(1,3), title = 'gtex scaffold w/o outlier') + theme(legend.position = "none")
ggsave('plots/gtex_scaffold_wo_outlier_dim13.png', width = 10, height = 8.72)

