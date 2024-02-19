library(tidyverse)
library(genefilter)
library(hgu133plus2.db)
library(e1071)
library(class)
library(GSVA)
library(genefilter)
library(affy)
library(AnnotationDbi)

# Load data
load(file = "../data/Breast cancer classification/CITBCMST.rdata") # CIT
load(file = "../data/Breast cancer classification/Bordet_paired_data/Bordet.rdata") # CIT
load("../data/Breast cancer classification/TCGA_RNA_seq/brca.Rdata") # PAM50
load("../data/Breast cancer classification/PAM50genes.Rdata") # PAM50

####################### CIT full ####################### 
idx <- findLargest(row.names(CIT_full), rowMedians(CIT_full), data="hgu133plus2")
matrix.exprs <- CIT_full[rownames(CIT_full) %in% idx,]
keys=keys(hgu133plus2.db)
symb=select(hgu133plus2.db, keys, "SYMBOL", "PROBEID" )
symb=symb[!duplicated(symb[,1]),]
row.names(matrix.exprs) <- symb[match(row.names(matrix.exprs), symb[,1]),][,2]
CIT_full_fl_genesymbol <- matrix.exprs[na.omit(row.names(matrix.exprs)),]

saveRDS(CIT_full_fl_genesymbol, '02_ss2bulk/rds/CIT_full_fl_genesymbol.rds')

####################### CIT core ####################### 
idx <- findLargest(row.names(CIT_core), rowMedians(CIT_core), data="hgu133plus2")
matrix.exprs <- CIT_core[rownames(CIT_core) %in% idx,]
keys=keys(hgu133plus2.db)
symb=select(hgu133plus2.db, keys, "SYMBOL", "PROBEID" )
symb=symb[!duplicated(symb[,1]),]
row.names(matrix.exprs) <- symb[match(row.names(matrix.exprs), symb[,1]),][,2]
CIT_core_fl_genesymbol <- matrix.exprs[na.omit(row.names(matrix.exprs)),]

saveRDS(CIT_core_fl_genesymbol, '02_ss2bulk/rds/CIT_core_fl_genesymbol.rds')

####################### PAM50 ####################### 

brca_pam50 # subtype/meta data
brca # expression data

