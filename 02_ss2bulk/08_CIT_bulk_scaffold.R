library(spaceRATScaffolds) 
library(spaceRAT)
library(tidyverse)
library(glue)
library(cowplot)

############## Versions of spaceRAT ##############
# remotes::install_github("shdam/spaceRAT", ref = 'groupScaling') # does not work 
# remotes::install_github("shdam/spaceRAT") # works
# packageDescription('spaceRAT')$GithubRef # Check what branch is installed 
##################################################

# Load data
CIT_full_fl_genesymbol <- readRDS('02_ss2bulk/rds/CIT_full_fl_genesymbol.rds')
CIT_core_fl_genesymbol <- readRDS('02_ss2bulk/rds/CIT_core_fl_genesymbol.rds')
load(file = "../data/Breast cancer classification/CITBCMST.rdata")

# Define CIT_classes as dataframe
CIT_classes_df <- as.data.frame(CIT_classes)

##############################################################
########################## CIT_FULL ########################## 
##############################################################

subset_deg <- TRUE
# pval_cutoff_list <- c(0.01, 0.05, 0.1) 
# lfc_cutoff_list <- c(1, 1.5, 2)
# n_genes_list <- c(1000, 500, 250)
# sort.by_list <- c('B', 'logFC', 'AveExpr', 't', 'P', 'none')

# for (pval_cutoff in pval_cutoff_list){
#   for (lfc_cutoff in lfc_cutoff_list){
#     for (n_genes in n_genes_list){
#       for (sort.by in sort.by_list){

# Best scaffold
pval_cutoff <- 0.05
lfc_cutoff <- 1.5
n_genes <- 250
sort.by <- 'B'

scaffold_full <- buildScaffold(object = CIT_full_fl_genesymbol, 
                               pheno = CIT_classes_df,
                               colname = "CIT_classes",
                               data = "exprs",
                               annotation = NULL,
                               pca_scale = TRUE,
                               subset_deg = subset_deg, 
                               pval_cutoff = pval_cutoff, 
                               lfc_cutoff = lfc_cutoff, 
                               n_genes = n_genes,
                               sort.by = sort.by
)

length(scaffold_full$DEgenes)

plotScaffold(scaffold_full, 
             title = glue("CIT full, subset_deg = {subset_deg}\n
                          pval: {pval_cutoff}, lfc: {lfc_cutoff}, Ngenes per group: {n_genes}, DEgenes: {length(scaffold_full$DEgenes)}, sort by: {sort.by}, "),
             dims = c(1,2)) + theme(legend.position = 'none')

ggsave(glue('02_ss2bulk/plots/08_CIT_full_scaffold_subset_deg_{subset_deg}_pval{pval_cutoff}_lfc{lfc_cutoff}_Ngenes{n_genes}_DEgenes{length(scaffold_full$DEgenes)}_sortby{sort.by}.png'))

#       }
#     }
#   }
# }


# Create an empty list to store plots

n_dims <- 4
P <- list()
count <- 0

for (dim1 in 1:n_dims) {
  
  for (dim2 in 1:n_dims) {
    
    if (dim1 == dim2 | dim1 > dim2){
      next
    }

    # Generate your plot using plotScaffold() function
    p <- plotScaffold(scaffold_full,
                      title = NULL,  # Set individual plot titles to NULL
                      dims = c(dim1, dim2)) + theme(legend.position = 'none')

    # Store the plot in list
    count <- count + 1
    P[[count]] <- p
  }
}

# Combine the plots into a grid
plot_grid(plotlist = P, ncol = n_dims-1, nrow = n_dims-2, labels = NULL)
ggsave(glue('02_ss2bulk/plots/08_CIT_full_scaffold_subset_deg_{subset_deg}_pval{pval_cutoff}_lfc{lfc_cutoff}_PC{n_dims}.png'),
       width = 13, height = 10, units = 'in'
       )



##############################################################
########################## CIT_core ########################## 
##############################################################

subset_deg = FALSE

scaffold_core <- buildScaffold(object = CIT_core_fl_genesymbol, 
                               pheno = CIT_classes_df,
                               colname = "CIT_classes",
                               data = "exprs",
                               annotation = NULL,
                               pca_scale = TRUE,
                               subset_deg = subset_deg)

plotScaffold(scaffold_core, 
             title = glue("CIT core, subset_deg = {subset_deg}"),
             dims = c(1,2)) + theme(legend.position = 'none')

ggsave(glue('02_ss2bulk/plots/08_CIT_core_scaffold_subset_deg_{subset_deg}.png'))

#####################################################################
########################### Best scaffold ########################### 
#####################################################################

subset_deg = TRUE
pval_cutoff <- 0.05
lfc_cutoff <- 1.5
n_genes <- 250
sort.by <- 'B'

best_scaffold <- buildScaffold(object = CIT_full_fl_genesymbol, 
                               pheno = CIT_classes_df,
                               colname = "CIT_classes",
                               data = "exprs",
                               annotation = NULL,
                               pca_scale = TRUE,
                               subset_deg = subset_deg, 
                               pval_cutoff = pval_cutoff, 
                               lfc_cutoff = lfc_cutoff, 
                               n_genes = n_genes,
                               sort.by = sort.by
)

best_scaffold_genes <- best_scaffold$DEgenes
core_genes <- rownames(CIT_core_fl_genesymbol)

table(core_genes %in% best_scaffold_genes)

saveRDS(best_scaffold, '02_ss2bulk/rds/CIT_best_scaffold.rds')


#####################################################################
##################### Best scaffold - remove samples ################ 
#####################################################################

samples_subset <- c('CIT_DSOA_034', 'CIT_DSOA_033', 'CIT_DSOA_035', 'CIT_DSOA_031', 'CIT_DSOA_535', 'CIT_DSOA_028')
# Expr
CIT_full_fl_genesymbol_remove_samples <- CIT_full_fl_genesymbol[, !colnames(CIT_full_fl_genesymbol) %in% samples_subset]
# Pheno
CIT_classes_df_remove_samples <- CIT_classes_df[!rownames(CIT_classes_df) %in% samples_subset ,]
CIT_classes_df_remove_samples <- as.data.frame(CIT_classes_df_remove_samples)
rownames(CIT_classes_df_remove_samples) <- colnames(CIT_full_fl_genesymbol_remove_samples)
colnames(CIT_classes_df_remove_samples) <- c('CIT_classes')
  
subset_deg = TRUE
pval_cutoff <- 0.05
lfc_cutoff <- 1.5
n_genes <- 250
sort.by <- 'B'

best_scaffold_remove_samples <- buildScaffold(object = CIT_full_fl_genesymbol_remove_samples, 
                                              pheno = CIT_classes_df_remove_samples,
                                              colname = "CIT_classes",
                                              data = "exprs",
                                              annotation = NULL,
                                              pca_scale = TRUE,
                                              subset_deg = subset_deg, 
                                              pval_cutoff = pval_cutoff, 
                                              lfc_cutoff = lfc_cutoff, 
                                              n_genes = n_genes,
                                              sort.by = sort.by
)

############# Project ############# 
# Pheno
project_samples_pheno <- CIT_classes_df[samples_subset ,]
project_samples_pheno <- as.data.frame(project_samples_pheno)
rownames(project_samples_pheno) <- samples_subset

# Expression
project_samples_expr <- CIT_full_fl_genesymbol[, samples_subset]

projectSample(scaffold = best_scaffold_remove_samples, 
              sample = project_samples_expr,
              pheno = project_samples_pheno,
              colname = 'project_samples_pheno',
              annotation = NULL,
              title = 'CIT samples projected to CIT scaffold')


