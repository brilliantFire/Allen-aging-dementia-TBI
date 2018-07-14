##
'
Feature Engineering & Data Preparation ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
13 Jul 2018 - Created

This script performs hierarchical clustering on expression levels for subsets of genes...

See http://aging.brain-map.org/overview/home for more details about the data.

R 3.4.1
'
##

'
Load libraries
'
library(edgeR)          # DGE library from Bioconductor
library(data.table)     # I/O
library(dplyr)          # Entering the TIDYVERSE!
library(VennDiagram)    # Venn diagrams with loads of graphical control
library(gplots)         # heatmaps, extracting from Venn diagrams
library(reshape2)       # split open and melt dataframes
library(ggplot2)        # other graphics

setwd('..')

'
Load data
'
fpkm_table <- readRDS(file='data/normalized_fpkm_matrix.Rds')
sample_info <- readRDS(file='data/sample_info.Rds')
brain_reg_sig_genes <- readRDS(file='data/brain_reg_sig_genes.Rds')
genes <- readRDS(file='data/genes.Rds')

hip_genes <- brain_reg_sig_genes$hip_genes
fwm_genes <- brain_reg_sig_genes$fwm_genes
pcx_genes <- brain_reg_sig_genes$pcx_genes
tcx_genes <- brain_reg_sig_genes$tcx_genes

sample_dementia_status <-  sample_info[, c('rnaseq_profile_id', 'act_demented')]

'
Heatmaps/clustering for filtered genes, all samples
'
# make a numeric matrix
fpkm_mat <- as.matrix(fpkm_table)
class(fpkm_mat) <- 'numeric'

# standardize
fpkm_standard_mat <- t(scale(t(fpkm_mat)))

# get all brain group sig genes
all_sig_genes <- unique(c(hip_genes, fwm_genes, pcx_genes, tcx_genes))

# heatmap for 30 random samples/columns
#png(filename='data/all_sig_brain_region_genes_heatmap_w_key_BIG.png', height = 1500, width = 1500)
set.seed(420)

# To specify a different liakage type: heatmap.2(...,hclustfun = function(x) hclust(x,method = 'centroid'),...)
heatmap.2(fpkm_standard_mat[all_sig_genes, sample(ncol(fpkm_standard_mat), 377)],
          hclustfun = function(x) hclust(x,method = 'complete'),
          distfun = function(x) dist(x,method='euclidean'),
          col=colorpanel(30, 
                         low = 'red', 
                         mid = 'white',
                         high = 'blue'),
          denscol='black',
          na.rm=TRUE,
          scale='none',
          trace = 'none', 
          dendrogram = 'both', 
          key = TRUE, 
          key.title = 'Key',
          cexRow = 0.4,
          cexCol = 0.4,
          margins = c(1.5,1.5),
          xlab='samples',
          ylab='genes',
          labRow=FALSE,
          labCol=FALSE)
#dev.off()

'
Exact tests for all genes, all samples
'
# Load & prep raw read counts
counts <- data.frame(read.csv('data/raw_read_counts.csv'))
rownames(counts) <- counts$gene_id
counts$gene_id <- NULL
colnames(counts) <- substring(colnames(counts), 2)
counts <- as.matrix(round(counts,0))
class(counts) <- 'numeric'

# get dementia status
counts_samples <- colnames(counts)
dementia_status <- sample_dementia_status$act_demented[sample_dementia_status$rnaseq_profile_id %in% counts_samples]

##
# DGEList() object
all_dgeList <- DGEList(counts, group = dementia_status)

# filter
all_filtered <- all_dgeList[rowSums(1e+06*all_dgeList$counts/expandAsMatrix(all_dgeList$samples$lib.size, dim(counts))>2)>=10,]
print(paste('Number of genes left after filtering: ', dim(all_filtered)[1]))

# calcNormFactors()
all_norm <- calcNormFactors(all_filtered)

# estimateGLMCommonDisp(), estimateGLMTrendedDisp(), estimateGLMTagwiseDisp()
all_norm <- estimateGLMCommonDisp(all_norm, verbose = TRUE)
all_norm <- estimateGLMTrendedDisp(all_norm)
all_norm <- estimateGLMTagwiseDisp(all_norm)

# Exact tests
all_test <- exactTest(all_norm, pair = c('No Dementia' , 'Dementia'), dispersion = 'auto')

# top tags
all_table <- topTags(all_test, n = 10000, p.value = 0.0001, adjust.method = 'BH', sort.by = 'PValue')
top_genes <- rownames(all_table$table)

# standardize
fpkm_standard <- t(scale(t(fpkm_table)))

# heatmap
heatmap.2(fpkm_standard[top_genes , sample(ncol(fpkm_standard), 377)],
          hclustfun = function(x) hclust(x, method = 'complete'),
          distfun = function(x) dist(x,method='euclidean'),
          col=colorpanel(20, 
                         low = 'red', 
                         mid = 'white',
                         high = 'blue'),
          na.rm=TRUE,
          scale='none',
          trace = 'none', 
          dendrogram = 'both', 
          key = TRUE, 
          key.title = 'Key',
          cexRow = 0.4,
          cexCol = 0.4,
          margins = c(1.5,1.5),
          xlab='samples',
          ylab='genes',
          labRow=FALSE,
          labCol=FALSE)

# gene list
sig_genes <- genes[genes$gene_entrez_id %in% top_genes,]
