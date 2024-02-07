library(spaceRATScaffolds)
library(spaceRAT)
library(tidyverse)

load(file="/home/people/s165827/data/Breast cancer classification/tcga_scaffold_v2.Rdata")
tcga_pheno <- as.data.frame(class_sample)
colnames(tcga_scaffold) <- paste0('sample-', 1:640)
rownames(tcga_pheno) <- colnames(tcga_scaffold)

scaffold <- getScaffold("GTEx.v1")
plotScaffold(scaffold)


################################################################################
################################## Baseline #################################### 
################################################################################

tcga_Scaffold_baseline <- buildScaffold(tcga_scaffold, 
                                        pheno = tcga_pheno,
                                        colname = 'class_sample',
                                        data = 'exprs', 
                                        subset_deg = TRUE,
                                        threshold = 10,
                                        pca_scale = TRUE)

plotScaffold(tcga_Scaffold_baseline, dims = c(1,2), title = 'tcga scaffold baseline') + theme(legend.position = "none")
ggsave('plots/tcga_scaffold_baseline_dim12.png', width = 10, height = 8.72)
plotScaffold(tcga_Scaffold_baseline, dims = c(2,3), title = 'tcga scaffold baseline') + theme(legend.position = "none")
ggsave('plots/tcga_scaffold_baseline_dim23.png', width = 10, height = 8.72)
plotScaffold(tcga_Scaffold_baseline, dims = c(1,3), title = 'tcga scaffold baseline') + theme(legend.position = "none")
ggsave('plots/tcga_scaffold_baseline_dim13.png', width = 10, height = 8.72)

################################################################################
################################### DE genes ################################### 
################################################################################

source('01_bluk2bluk/functions.R')

DEsubset <- c('rds/tcga_DE_L2FC_0.3_pval_0.05.rds',
              'rds/tcga_DE_pval_0.05.rds',
              'rds/tcga_20DEgenes.rds')

lapply(DEsubset, 
       subset_scaffold, 
       basline_scaffold = tcga_scaffold, 
       class_sample = class_sample, 
       extra_name = "_TCGA")

# subset_scaffold(rds_path = 'rds/tcga_DE_L2FC_0.3_pval_0.05.rds',
#                 basline_scaffold = tcga_scaffold, 
#                 class_sample = class_sample, 
#                 extra_name = "_TCGA")



