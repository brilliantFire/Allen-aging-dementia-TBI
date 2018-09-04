##
'
Feature Engineering & Data Preparation ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
13 Jul 2018 - Created
31 Aug 2018 - clValid

This script performs hierarchical clustering on expression levels for subsets of genes...

See http://aging.brain-map.org/overview/home for more details about the data.

R 3.4.3
'
##

'
Load libraries
'
library(edgeR)          # DGE library from Bioconductor
library(data.table)     # I/O
library(dplyr)          # Entering the TIDYVERSE!
library(gplots)         # heatmaps, extracting from Venn diagrams
library(reshape2)       # split open and melt dataframes
library(cluster)

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


set.seed(420)

# To specify a different linkage type: heatmap.2(...,hclustfun = function(x) hclust(x,method = 'centroid'),...)

# heatmap for TCx sig genes
tcx_hm <- heatmap.2(fpkm_standard_mat[tcx_genes, sample(ncol(fpkm_standard_mat), 377)],
                    hclustfun = function(x) hclust(x,method = 'ward.D2'),
                    distfun = function(x) as.dist(1-cor(t(x))), # 1 minus the Pearson corr distance
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

# pull gene dendrogram
tcx_dendro <- tcx_hm$rowDendrogram

# brute force: distance matrix
tcx_dist <- as.dist(1-cor(t(fpkm_standard_mat[tcx_genes, sample(ncol(fpkm_standard_mat), 377)])))
# cluster
tcx_cluster <- hclust(tcx_dist,method = 'complete')

'
Clustering Optimization
'
source('https://bioconductor.org/biocLite.R')
biocLite('GenomicAlignments')

# 'org.Hs.eg.db' - Entrez Gene identifier-based annotation
# 'mygene' - for constructing lists of functional classes
biocLite(c('Biobase', 'annotate', 'GO.db', 'org.Hs.eg.db', 'mygene'))

# install clValid from CRAN
install.packages('clValid', repo='https://CRAN.R-project.org/')

# load libraries
library(Biobase)
library(annotate)
library(GO.db)
library(org.Hs.eg.db)
library(mygene)
library(clValid)

# HIP FPKM subset
hip_samples <- sample_info$rnaseq_profile_id[which(sample_info$structure_acronym == 'HIP')]
hip_data <- fpkm_standard_mat[hip_genes, colnames(fpkm_standard_mat) %in% hip_samples]

'
Gene ontology lists from mygene.info database
'
# gene ontology info for Entrez IDs
hip_res <- queryMany(rownames(hip_data), scopes='entrezgene', fields=c('go', 'symbol'), species='human')

# take a look at the first few entries
head(hip_res, 10)

# pull out Entrez Gene IDs & symbols
hip_entrez <- hip_res$query
hip_symbols <- hip_res$symbol

# take a look at a few individual gene ontology records; CC for gene #1 (Entrez Gene #25937, WWTR1)
hip_res[1, 'go.CC'][[1]]
# MF for gene #1
hip_res[1, 'go.MF'][[1]]
# BP for gene #1
hip_res[1, 'go.BP'][[1]]

# for gene #70 (Entrez Gene #23580, CDC42EP4)
hip_res[70, 'go.CC'][[1]]
hip_res[70, 'go.MF'][[1]]
hip_res[70, 'go.BP'][[1]]

# extracts most common functional category for each gene for three ontologies:
#    CC = cellular component
#    MF = molecular function
#    BP = biological process
hip_cc <- NULL
hip_mf <- NULL
hip_bp <- NULL
for(i in 1:length(hip_genes)){
    # each if-else checks each 'term' list's length and, if zero, marks as 'unknown'.
    # if not zero, returns the most common function
    if(length(hip_res[i, 'go.CC'][[1]]$term) == 0){
        a_hip_cc <- 'unknown'
    }
    else{
        a_hip_cc <- tail(names(sort(table(hip_res[i, 'go.CC'][[1]]$term))), 1)
    } 
    hip_cc <- c(hip_cc, a_hip_cc)
    
    if(length(hip_res[i, 'go.MF'][[1]]$term) == 0){
        a_hip_mf <- 'unknown'
    }
    else{
        a_hip_mf <- tail(names(sort(table(hip_res[i, 'go.MF'][[1]]$term))), 1)
    } 
    hip_mf <- c(hip_mf, a_hip_mf)
        
    if(length(hip_res[i, 'go.BP'][[1]]$term) == 0){
        a_hip_bp <- 'unknown'
    }
    else{
        a_hip_bp <- tail(names(sort(table(hip_res[i, 'go.BP'][[1]]$term))), 1)
    } 
    hip_bp <- c(hip_bp, a_hip_bp)
}

hip_fc <- data.frame(cbind(hip_cc, hip_mf, hip_bp))
rownames(hip_fc) <- hip_entrez
colnames(hip_fc) <- c('CC', 'MF', 'BP')

hip_fc_cc <- tapply(rownames(hip_data), hip_fc$CC, c)

'
Clustering Optimization
'
hip_valid <- clValid(hip_data, 2:50, 
                     clMethods = c('hierarchical', 'kmeans', 'diana', 'pam', 'clara', 'model'), 
                     validation = c('internal', 'stability', 'biological'),
                     maxitems = length(hip_genes),
                     metric = 'euclidean',
                     method = 'complete',
                     annotation = hip_fc_cc,
                     verbose = TRUE)
summary(hip_valid)
optimalScores(hip_valid)

##########
plot(test, measure = "Dunn", legendLoc = "topleft")

hc <- clusters(test, "hierarchical")
fc_labels <- factor(hip_fc$CC)
plot(hc, labels = fc_labels)



seven <- cutree(hc, 7)
hip_clusters <- xtabs(~hip_fc$CC + seven)




















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
