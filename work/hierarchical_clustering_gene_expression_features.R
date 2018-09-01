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

# install clValid package from CRAN repo
install.packages('clValid', repos='http://cran.us.r-project.org')
install.packages('kohonen', repos='http://cran.us.r-project.org')
install.packages('mclust', repos='http://cran.us.r-project.org')
library(clValid)
library(kohonen)
library(mclust)

# get annotation packages from Bioconductor
source('http://bioconductor.org/biocLite.R')
biocLite(c('Biobase', 'annotate', 'GO.db'))
biocLite('org.Hs.eg.db') # Entrez Gene identifier-based annotation
biocLite('mygene')       # for functional classes

# HIP FPKM subset
hip_samples <- sample_info$rnaseq_profile_id[which(sample_info$structure_acronym == 'HIP')]
hip_data <- fpkm_standard_mat[hip_genes, colnames(fpkm_standard_mat) %in% hip_samples]

# trying out clValid
library(Biobase)
library(annotate)
library(GO.db)
library(org.Hs.eg.db)
library(mygene)

# gene ontology info for Entrez IDs
hip_res <- queryMany(rownames(hip_data), scopes='entrezgene', fields=c('entrezgene', 'go', 'symbol'), species='human')

hip_res$go.CC[which(hip_res$go.CC == 'NULL')] <- 'term'

hip_res$go.CC <- lapply(hip_res$go.CC, function(z) { z[ lengths(z) == 0 ] <- list(term = 'unknown'); z; })


# extract most common FC for a single gene
tail(names(sort(table(res[100, 'go.MF'][[1]]$term))), 1)

hip_cc <- NULL
hip_mf <- NULL
hip_bp <- NULL
for(i in 1:length(hip_genes)){
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
rownames(hip_fc) <- hip_res$entrezgene
colnames(hip_fc) <- c('CC', 'MF', 'BP')

hip_fc$CC[hip_fc$CC == 'cytosol'] = 'cytoplasm'
fc <- tapply(rownames(hip_data), hip_fc$CC, c)
fc <- fc[!names(fc) %in% c('unknown')]


test <- clValid(hip_data, 2:20, 
                clMethods = c('hierarchical'), 
                validation = 'biological',
                maxitems = length(hip_genes),
                metric = 'correlation',
                method = 'complete',
                annotation = fc,
                #GOcategory = 'all',
                verbose = TRUE)
summary(test)

plot(test, measure = "BHI", legendLoc = "topleft")

hc <- clusters(test, "hierarchical")
fc_labels <- factor(hip_fc$CC)
plot(hc, labels = fc_labels)

optimalScores(test)

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
