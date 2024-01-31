
library(spaceRATScaffolds)
library(spaceRAT)
library(tidyverse)

listScaffolds()

gtex <- getScaffold('GTEx.v1')
plotScaffold(gtex, plot_mode = 'tiny_label')
gtex$pca$center

tcga_testis <- read.csv('../../Speciale/tcga_testis', sep = '\t' )
rownames(tcga_testis) <- NULL
tcga_testis <- column_to_rownames(tcga_testis, var = 'sample')
head(tcga_testis)

projectSample(scaffold = gtex,
              sample = tcga_testis, 
              title = 'Testis (non-normalized) samples projected onto scaffold PCA.')

ggsave('../../Speciale/out/tcga_testis.png',
       width = 16,
       height = 10)

tcga_testis_normalized <- read.csv('../../Speciale/tcga_testis_normalized', sep = '\t' )
rownames(tcga_testis_normalized) <- NULL
tcga_testis_normalized <- column_to_rownames(tcga_testis_normalized, var = 'sample')
head(tcga_testis_normalized)

projectSample(scaffold = gtex,
              sample = tcga_testis_normalized, 
              title = 'Testis (normalized) samples projected onto scaffold PCA.')

ggsave('../../Speciale/out/tcga_testis_normalized.png',
       width = 16,
       height = 10)

tcga_testis_fpkm <- read.csv('../../Speciale/tcga_testis_fpkm.tsv', sep = '\t' )
rownames(tcga_testis_fpkm) <- NULL
tcga_testis <- column_to_rownames(tcga_testis_fpkm, var = 'Ensembl_ID')
head(tcga_testis_fpkm)
projectSample(scaffold = gtex,
              sample = tcga_testis_fpkm, 
              title = 'Testis (FPKM) samples projected onto scaffold PCA.')
ggsave('../../Speciale/out/tcga_testis_fpkm.png',
       width = 16,
       height = 10)


tcga_DLBC <- read.csv('../../Speciale/tcga_DLBC', sep = '\t' )
rownames(tcga_DLBC) <- NULL
tcga_DLBC <- column_to_rownames(tcga_DLBC, var = 'sample')
head(tcga_DLBC)

projectSample(scaffold = gtex,
              sample = tcga_DLBC)

tcga_DLBC_normalized <- read.csv('../../Speciale/tcga_DLBC_normalized', sep = '\t' )
rownames(tcga_DLBC_normalized) <- NULL
tcga_DLBC_normalized <- column_to_rownames(tcga_DLBC_normalized, var = 'sample')
head(tcga_DLBC_normalized)

projectSample(scaffold = gtex,
              sample = tcga_DLBC_normalized)


