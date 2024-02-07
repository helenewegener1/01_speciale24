library(spaceRATScaffolds)
library(spaceRAT)
library(tidyverse)

load(file="/home/people/s165827/data/Breast cancer classification/gtex.Rdata")
source('01_bluk2bluk/functions.R')

################################################################################
############################### Subset samples ################################# 
################################################################################

gtex_top10_samples <- subset_sample(scaffold = gtex_scaffold, class_sample = class_sample, n_samples = 10, rank = FALSE)
gtex_top10_samples_scaffold <- gtex_top10_samples[[1]]
gtex_top10_samples_pheno <- gtex_top10_samples[[2]]

gtex_top10_samples_Scaffold <- buildScaffold(gtex_top10_samples_scaffold, 
                                             pheno = gtex_top10_samples_pheno,
                                             colname = 'class_sample',
                                             data = 'exprs', 
                                             subset_deg = FALSE,
                                             threshold = 10,
                                             # add_umap = TRUE,
                                             pca_scale = TRUE)

p1 <- plotScaffold(gtex_top10_samples_Scaffold, dims = c(1,2), title = 'gtex scaffold top 10 sample') + theme(legend.position = "none")
ggsave('plots/gtex_scaffold_top10sample_dim12.png', width = 10, height = 8.72)
p2 <- plotScaffold(gtex_top10_samples_Scaffold, dims = c(2,3), title = 'gtex scaffold top 10 sample') + theme(legend.position = "none")
ggsave('plots/gtex_scaffold_top10sample_dim23.png', width = 10, height = 8.72)
p3 <- plotScaffold(gtex_top10_samples_Scaffold, dims = c(1,3), title = 'gtex scaffold top 10 sample') + theme(legend.position = "none")
ggsave('plots/gtex_scaffold_top10sample_dim13.png', width = 10, height = 8.72)

return_plot <- wrap_plots(p1, p2, p3, ncol = 2)
ggsave(glue('plots/gtex_scaffold_top10sample_wrapped.png'), width = 15, height = 8)

################################################################################
###################### Subset samples, subset_deg = TRUE ####################### 
################################################################################

gtex_top10_samples_subset_deg_Scaffold <- buildScaffold(gtex_top10_samples_scaffold, 
                                                        pheno = gtex_top10_samples_pheno,
                                                        colname = 'class_sample',
                                                        data = 'exprs', 
                                                        subset_deg = TRUE,
                                                        threshold = 10,
                                                        # add_umap = TRUE,
                                                        pca_scale = TRUE)

plotScaffold(gtex_top10_samples_subset_deg_Scaffold, dims = c(1,2), title = 'gtex scaffold top 10 sample subset_deg = TRUE') + theme(legend.position = "none")
ggsave('plots/gtex_scaffold_top10sample_subset_deg_dim12.png', width = 10, height = 8.72)
plotScaffold(gtex_top10_samples_subset_deg_Scaffold, dims = c(2,3), title = 'gtex scaffold top 10 sample subset_deg = TRUE') + theme(legend.position = "none")
ggsave('plots/gtex_scaffold_top10sample_subset_deg_dim23.png', width = 10, height = 8.72)
plotScaffold(gtex_top10_samples_subset_deg_Scaffold, dims = c(1,3), title = 'gtex scaffold top 10 sample subset_deg = TRUE') + theme(legend.position = "none")
ggsave('plots/gtex_scaffold_top10sample_subset_deg_dim13.png', width = 10, height = 8.72)


################################################################################
############################ Subset samples - rank ############################# 
################################################################################

# gtex_top10_samples_rank <- subset_sample(scaffold = gtex_scaffold, class_sample = class_sample, n_samples = 10, rank = TRUE)
# gtex_top10_samples_rank_scaffold <- gtex_top10_samples_rank[[1]]
# gtex_top10_samples_rank_pheno <- gtex_top10_samples_rank[[2]]

################################################################################
################################### DE genes ################################### 
################################################################################

DE_subsets <- c('rds/DE_L2FC_0.3_pval_0.05.rds', 
                'rds/DE_L2FC_0.2_pval_0.05.rds', 
                'rds/DE_L2FC_0.1_pval_0.05.rds', 
                'rds/DE_pval_0.05.rds', 
                'rds/20DEgenes.rds')

lapply(DE_subsets, 
       subset_scaffold, 
       basline_scaffold = gtex_scaffold, 
       class_sample = class_sample)

################################################################################
############################ DE genes - w/o outlier ############################ 
################################################################################

readRDS('rds/gtex_scaffold_wo_outlier.rds')
readRDS('rds/gtex_pheno_wo_outlier.rds')

class_sample_wo_outlier <- gtex_pheno_wo_outlier$class_sample %>% unlist() 

lapply(DE_subsets, 
       subset_scaffold, 
       basline_scaffold = gtex_scaffold_wo_outlier, 
       class_sample = class_sample_wo_outlier, 
       extra_name = "_wo_outlier")


################################################################################
########################## DE genes + top 10 samples ########################### 
################################################################################

class_sample_top10samples <- gtex_top10_samples_pheno$class_sample %>% unlist()

lapply(subsets, 
       subset_scaffold, 
       basline_scaffold = gtex_top10_samples_scaffold, 
       class_sample = class_sample_top10samples, 
       extra_name = "_top10samples")














################################################################################
######################### Subset samples + DE genes ############################
################################################################################

gtex_20DEgenes_top10_samples <- subset_sample(scaffold = gtex_scaffold_20DEgenes_subset, class_sample = class_sample, n_samples = 10, rank = FALSE)
gtex_20DEgenes_top10_samples_scaffold <- gtex_20DEgenes_top10_samples[[1]]
gtex_20DEgenes_top10_samples_pheno <- gtex_20DEgenes_top10_samples[[2]]

gtex_20DEgenes_top10_samples_Scaffold <- buildScaffold(gtex_20DEgenes_top10_samples_scaffold, 
                                                       pheno = gtex_20DEgenes_top10_samples_pheno,
                                                       colname = 'class_sample',
                                                       data = 'exprs', 
                                                       subset_deg = FALSE,
                                                       threshold = 10,
                                                       # add_umap = TRUE,
                                                       pca_scale = TRUE)

plotScaffold(gtex_20DEgenes_top10_samples_Scaffold, dims = c(1,2), title = 'gtex scaffold top 10 sample + top 20 DE genes') + theme(legend.position = "none")
ggsave('plots/gtex_scaffold_DEsubset_top10sample_dim12.png', width = 10, height = 8.72)
plotScaffold(gtex_20DEgenes_top10_samples_Scaffold, dims = c(2,3), title = 'gtex scaffold top 10 sample + top 20 DE genes') + theme(legend.position = "none")
ggsave('plots/gtex_scaffold_DEsubset_top10sample_dim23.png', width = 10, height = 8.72)
plotScaffold(gtex_20DEgenes_top10_samples_Scaffold, dims = c(1,3), title = 'gtex scaffold top 10 sample + top 20 DE genes') + theme(legend.position = "none")
ggsave('plots/gtex_scaffold_DEsubset_top10sample_dim13.png', width = 10, height = 8.72)


