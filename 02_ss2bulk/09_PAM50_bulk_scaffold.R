library(spaceRATScaffolds) 
library(spaceRAT)
library(tidyverse)
library(glue)
library(cowplot)

# Load data
load("../data/Breast cancer classification/TCGA_RNA_seq/brca.Rdata") # PAM50
load("../data/Breast cancer classification/PAM50genes.Rdata") # PAM50

# Wrangle pheno data 
brca_pam50 <- brca_pam50[match(colnames(brca), brca_pam50$sample),] # keep onyl sample that is also in expression data
brca_pam50_new <- as.data.frame(brca_pam50[brca_pam50$sample %in% colnames(brca),]) 
rownames(brca_pam50_new) <- NULL
brca_pam50_new <- brca_pam50_new %>% column_to_rownames('sample')

##############################################################
######################### PAM50 FULL ######################### 
##############################################################

subset_deg <- TRUE
# pval_cutoff_list <- c(0.01, 0.05, 0.1) 
pval_cutoff_list <- c(0.05)
lfc_cutoff_list <- c(2.5)
n_genes_list <- c(50, 20)
# sort.by_list <- c('B', 'logFC', 'AveExpr', 't', 'P', 'none')
sort.by_list <- c('B')

for (pval_cutoff in pval_cutoff_list){
  for (lfc_cutoff in lfc_cutoff_list){
    for (n_genes in n_genes_list){
      for (sort.by in sort.by_list){
        
        scaffold_full <- buildScaffold(object = brca, 
                                       pheno = brca_pam50_new,
                                       colname = "PAM50",
                                       data = "exprs",
                                       annotation = NULL,
                                       pca_scale = TRUE,
                                       subset_deg = subset_deg, 
                                       pval_cutoff = pval_cutoff, 
                                       lfc_cutoff = lfc_cutoff, 
                                       n_genes = n_genes,
                                       sort.by = sort.by
        )
        
        # length(scaffold_full$DEgenes)
        
        plotScaffold(scaffold_full, 
                     title = glue("PAM50 full, subset_deg = {subset_deg}\n
                                  pval: {pval_cutoff}, lfc: {lfc_cutoff}, Ngenes per group: {n_genes}, DEgenes: {length(scaffold_full$DEgenes)}, sort by: {sort.by}, "),
                     dims = c(1,2)) + theme(legend.position = 'none')
        
        ggsave(glue('02_ss2bulk/plots/09_PAM50_full_scaffold_subset_deg_{subset_deg}_pval{pval_cutoff}_lfc{lfc_cutoff}_Ngenes{n_genes}_DEgenes{length(scaffold_full$DEgenes)}_sortby{sort.by}.png'))
        
      }
    }
  }
}

##############################################################
############################ PAM50 ########################### 
##############################################################

subset_deg = TRUE

# Subset full expression set to 
pam50_data <- as.matrix(brca[pam50_genes, ])

scaffold_core <- buildScaffold(object = pam50_data, 
                               pheno = brca_pam50_new,
                               colname = "PAM50",
                               data = "exprs",
                               annotation = NULL,
                               pca_scale = TRUE,
                               subset_deg = subset_deg,
                               pval_cutoff = pval_cutoff,
                               lfc_cutoff = lfc_cutoff
)

plotScaffold(scaffold_core, 
             title = glue("PAM50 core, subset_deg = {subset_deg} \n
                              pval: {pval_cutoff}, lfc: {lfc_cutoff}"),
             dims = c(1,2)) + theme(legend.position = 'none')

ggsave(glue('02_ss2bulk/plots/09_PAM50_core_scaffold_subset_deg_{subset_deg}_pval{pval_cutoff}_lfc{lfc_cutoff}.png'))


# plotScaffold(scaffold_core, 
#              title = glue("PAM50 core, subset_deg = {subset_deg}"),
#              dims = c(1,2)) + theme(legend.position = 'none')
# 
# ggsave(glue('02_ss2bulk/plots/09_PAM50_core_scaffold_subset_deg_{subset_deg}.png'))


##############################################################
######################## Overlap genes ####################### 
##############################################################

best_scaffold <- buildScaffold(object = brca, 
                               pheno = brca_pam50_new,
                               colname = "PAM50",
                               data = "exprs",
                               annotation = NULL,
                               pca_scale = TRUE,
                               subset_deg = TRUE, 
                               pval_cutoff = 0.05, 
                               lfc_cutoff = 1.5, 
                               n_genes = 20,
                               sort.by = 'B'
)

best_scaffold_genes <- best_scaffold$DEgenes

table(pam50_genes %in% best_scaffold_genes)

saveRDS(best_scaffold, '02_ss2bulk/rds/PAM50_best_scaffold.rds')

# plotScaffold(best_scaffold, 
#              title = glue("PAM50 Scaffold, pval_cutoff = 0.05, lfc_cutoff = 1.5, n_genes = 20, sort.by = B"),
#              dims = c(2,3)) + theme(legend.position = 'none')
# 
# ggsave('02_ss2bulk/plots/09_PAM50_best_scaffold_dim23.png')
 
