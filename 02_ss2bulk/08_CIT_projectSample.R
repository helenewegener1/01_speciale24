library(spaceRATScaffolds) 
library(spaceRAT)
library(tidyverse)
library(glue)
library(cowplot)

# Load data
best_scaffold <- readRDS('02_ss2bulk/rds/CIT_best_scaffold.rds')
load(file = "../data/Breast cancer classification/CITBCMST.rdata")
CIT_full_fl_genesymbol <- readRDS('02_ss2bulk/rds/CIT_full_fl_genesymbol.rds')

# Define CIT_classes as dataframe
CIT_classes_df <- as.data.frame(CIT_classes)

##############################################################
###################### Sample to project ##################### 
##############################################################

# Pheno
samples_subset <- c('CIT_DSOA_034', 'CIT_DSOA_033', 'CIT_DSOA_035', 'CIT_DSOA_031', 'CIT_DSOA_535', 'CIT_DSOA_028')
project_samples_pheno <- CIT_classes_df[samples_subset ,]
project_samples_pheno <- as.data.frame(project_samples_pheno)
rownames(project_samples_pheno) <- samples_subset

# Expression
project_samples_expr <- CIT_full_fl_genesymbol[, samples_subset]

saveRDS(project_samples_expr, '02_ss2bulk/rds/CIT_project_samples_expr.rds')
saveRDS(project_samples_pheno, '02_ss2bulk/rds/CIT_project_samples_pheno.rds')

##############################################################
######################## Project CIT ######################### 
##############################################################

projectSample(scaffold = best_scaffold, 
              sample = project_samples_expr,
              pheno = project_samples_pheno,
              colname = 'project_samples_pheno',
              title = 'CIT samples projected to CIT scaffold')

ggsave('02_ss2bulk/plots/08_CIT_project_CIT_samples_to_CIT_scaffold.png')

##############################################################
####################### Project PAM50 ######################## 
##############################################################

PAM50_project_samples_expr <- readRDS('02_ss2bulk/rds/PAM50_project_samples_expr.rds')
PAM50_project_samples_pheno <- readRDS('02_ss2bulk/rds/PAM50_project_samples_pheno.rds')

projectSample(scaffold = best_scaffold, 
              sample = PAM50_project_samples_expr,
              pheno = PAM50_project_samples_pheno,
              colname = 'project_samples_pheno',
              title = 'PAM50 samples projected to CIT scaffold')

ggsave('02_ss2bulk/plots/08_CIT_project_PAM50_samples_to_CIT_scaffold.png')
