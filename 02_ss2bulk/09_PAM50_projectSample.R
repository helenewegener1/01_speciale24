library(spaceRATScaffolds) 
library(spaceRAT)
library(tidyverse)
library(glue)
library(cowplot)

# Load data
load("../data/Breast cancer classification/TCGA_RNA_seq/brca.Rdata") # PAM50
# load("../data/Breast cancer classification/PAM50genes.Rdata") # PAM50

best_scaffold <- readRDS('02_ss2bulk/rds/PAM50_best_scaffold.rds')

# Wrangle pheno data
brca_pam50 <- brca_pam50[match(colnames(brca), brca_pam50$sample),] # keep onyl sample that is also in expression data
brca_pam50_new <- as.data.frame(brca_pam50[brca_pam50$sample %in% colnames(brca),])
rownames(brca_pam50_new) <- NULL
brca_pam50_new <- brca_pam50_new %>% column_to_rownames('sample')

##############################################################
###################### Sample to project ##################### 
##############################################################

# Pheno
samples_subset <- c('TCGA-C8-A1HL-01', 'TCGA-A2-A3XX-01', 'TCGA-AC-A3QP-01', 'TCGA-A8-A0A7-01', 'TCGA-LL-A441-01')
project_samples_pheno <- brca_pam50_new[samples_subset ,]
project_samples_pheno <- as.data.frame(project_samples_pheno)
rownames(project_samples_pheno) <- samples_subset

# Expression
project_samples_expr <- brca[, samples_subset]

saveRDS(project_samples_expr, '02_ss2bulk/rds/PAM50_project_samples_expr.rds')
saveRDS(project_samples_pheno, '02_ss2bulk/rds/PAM50_project_samples_pheno.rds')
  
##############################################################
######################## Project PAM50 ####################### 
##############################################################

projectSample(scaffold = best_scaffold, 
              sample = project_samples_expr,
              pheno = project_samples_pheno,
              colname = 'project_samples_pheno', 
              annotation = NULL,
              title = 'PAM50 samples projected to PAM50 scaffold')

ggsave('02_ss2bulk/plots/09_PAM50_project_PAM50_samples_to_PAM50_scaffold.png')

##############################################################
######################## Project CIT ######################### 
##############################################################

CIT_project_samples_expr <- readRDS('02_ss2bulk/rds/CIT_project_samples_expr.rds')
CIT_project_samples_pheno <- readRDS('02_ss2bulk/rds/CIT_project_samples_pheno.rds')

projectSample(scaffold = best_scaffold, 
              sample = CIT_project_samples_expr,
              pheno = CIT_project_samples_pheno,
              colname = 'project_samples_pheno', 
              annotation = NULL,
              title = 'CIT samples projected to PAM50 scaffold')

ggsave('02_ss2bulk/plots/09_PAM50_project_CIT_samples_to_PAM50_scaffold.png')






