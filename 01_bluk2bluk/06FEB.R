# Update spaceRAT 06FEB

library(spaceRATScaffolds)
library(spaceRAT)
library(tidyverse)

load(file="/home/people/s165827/data/Breast cancer classification/gtex.Rdata")
load(file='/home/people/s165827/data/Breast cancer classification/Bordet_paired_data/Bordet_RNA_gene_tpm.Rdata')
gtex_pheno <- as.data.frame(class_sample)
rownames(gtex_pheno) <- colnames(gtex_scaffold)

source('01_bluk2bluk/functions.R')
gtex_top10_samples <- subset_sample(scaffold = gtex_scaffold, class_sample = class_sample, n_samples = 10, rank = FALSE)
gtex_top10_samples_scaffold <- gtex_top10_samples[[1]]
gtex_top10_samples_pheno <- gtex_top10_samples[[2]]

################################################################################
################################## Baseline ################################### 
################################################################################

# for (n_genes in c(1200, 1700, 2000)){
  
n_genes = 1000
gtex_Scaffold_baseline <- buildScaffold(gtex_scaffold, 
                                        pheno = gtex_pheno,
                                        colname = 'class_sample',
                                        data = 'exprs', 
                                        subset_deg = TRUE,
                                        threshold = 10,
                                        pca_scale = TRUE,
                                        # add_umap = TRUE,
                                        n_genes = n_genes)

plotScaffold(gtex_Scaffold_baseline, dims = c(1,2), title = glue('gtex scaffold baseline {n_genes}')) + theme(legend.position = "none")
ggsave(glue('plots_FEB06/{n_genes}genes_baseline_dim12.png'), width = 10, height = 8.72)
plotScaffold(gtex_Scaffold_baseline, dims = c(2,3), title = glue('gtex scaffold baseline {n_genes}')) + theme(legend.position = "none")
ggsave(glue('plots_FEB06/{n_genes}genes_baseline_dim23.png'), width = 10, height = 8.72)
plotScaffold(gtex_Scaffold_baseline, dims = c(1,3), title = glue('gtex scaffold baseline {n_genes}')) + theme(legend.position = "none")
ggsave(glue('plots_FEB06/{n_genes}genes_baseline_dim13.png'), width = 10, height = 8.72)

####### Projecting samples in ####### 
# projectSample(scaffold = gtex_Scaffold_baseline,
#               sample = Bordet_RNA_gene_tpm, 
#               title = 'Bordet samples projected. N genes: {n_genes}', 
#               dims = c(2,3)) + theme(legend.position = "none")
# ggsave(glue('plots_FEB06/Breast_projected_{n_genes}genes_baseline_dim23.png'), width = 10, height = 8.72)

n_genes = 1000
gtex_top10_samples_Scaffold <- buildScaffold(gtex_top10_samples_scaffold, 
                                             pheno = gtex_top10_samples_pheno,
                                             colname = 'class_sample',
                                             data = 'exprs', 
                                             subset_deg = TRUE,
                                             threshold = 10,
                                             add_umap = TRUE,
                                             pca_scale = TRUE,
                                             n_genes = n_genes)

plotScaffold(gtex_top10_samples_Scaffold, dims = c(1,2), title = glue('gtex scaffold top10samples {n_genes}')) + theme(legend.position = "none")
ggsave(glue('plots_FEB06/{n_genes}genes_top10sample_dim12.png'), width = 10, height = 8.72)
plotScaffold(gtex_top10_samples_Scaffold, dims = c(2,3), title = glue('gtex scaffold top10samples {n_genes}')) + theme(legend.position = "none")
ggsave(glue('plots_FEB06/{n_genes}genes_top10sample_dim23.png'), width = 10, height = 8.72)
plotScaffold(gtex_top10_samples_Scaffold, dims = c(1,3), title = glue('gtex scaffold top10samples {n_genes}')) + theme(legend.position = "none")
ggsave(glue('plots_FEB06/{n_genes}genes_top10sample_dim13.png'), width = 10, height = 8.72)

####### Projecting samples in ####### 
projectSample(scaffold = gtex_top10_samples_Scaffold,
              sample = Bordet_RNA_gene_tpm, 
              title = glue('Bordet samples projected. N genes: {n_genes}'), 
              dims = c(2,3)) + theme(legend.position = "none")
ggsave(glue('plots_FEB06/Breast_projected_{n_genes}genes_top10sample_dim23.png'), width = 10, height = 8.72)

# }




# DE_subsets <- c('rds/DE_L2FC_0.3_pval_0.05.rds', 
#                 'rds/DE_pval_0.05.rds', 
#                 'rds/20DEgenes.rds')
# 
# lapply(DE_subsets, 
#        subset_scaffold, 
#        basline_scaffold = gtex_scaffold, 
#        class_sample = class_sample,
#        extra_name = '_NEW')

######



gtex_top10_samples_Scaffold <- buildScaffold(gtex_top10_samples_scaffold, 
                                             pheno = gtex_top10_samples_pheno,
                                             colname = 'class_sample',
                                             data = 'exprs', 
                                             subset_deg = TRUE,
                                             threshold = 10,
                                             add_umap = TRUE,
                                             pca_scale = TRUE,
                                             n_genes = 200)

plotScaffold(gtex_top10_samples_Scaffold, dims = c(1,2), title = 'gtex scaffold top 10 sample') + theme(legend.position = "none")
ggsave('plots/FEB06_200genes_gtex_scaffold_top10sample_dim12.png', width = 10, height = 8.72)
# plotScaffold(gtex_top10_samples_Scaffold, dims = c(2,3), title = 'gtex scaffold top 10 sample') + theme(legend.position = "none")
# ggsave('NEW_plots/gtex_scaffold_top10sample_dim23.png', width = 10, height = 8.72)
# plotScaffold(gtex_top10_samples_Scaffold, dims = c(1,3), title = 'gtex scaffold top 10 sample') + theme(legend.position = "none")
# ggsave('NEW_plots/gtex_scaffold_top10sample_dim13.png', width = 10, height = 8.72)

#####
class_sample_top10samples <- gtex_top10_samples_pheno$class_sample %>% unlist()

DE_subsets <- c('rds/DE_L2FC_0.3_pval_0.05.rds', 
                'rds/DE_L2FC_0.2_pval_0.05.rds', 
                'rds/DE_L2FC_0.1_pval_0.05.rds', 
                'rds/DE_pval_0.05.rds', 
                'rds/20DEgenes.rds')

lapply(DE_subsets, 
       subset_scaffold, 
       basline_scaffold = gtex_top10_samples_scaffold, 
       class_sample = class_sample_top10samples, 
       extra_name = "_top10samples_NEW")





