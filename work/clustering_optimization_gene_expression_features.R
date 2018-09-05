##
'
Feature Engineering & Data Preparation ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
13 Jul 2018 - Created
31 Aug 2018 - clValid tests
02 Sep 2018 - mygene annotation list function

This script uses the clValid package to run clustering experiments on normalized
FPKM gene expression values for each of the four regions of the brain in the 
dataset. It uses `mygene` to access gene ontology info in the `mygene.info`
ElasticSearch database and creates three ontology lists for each brain region:
    1. Cellular Component (CC)
    2. Molecular Function (MF)
    3. Biological Process (BP)
Clustering algorithms tested include agglomerative hierarchical clustering, CLARA,
and mixture model-based clustering. The "best" clustering strategies are initially
chosen on the basis of "internal" metrics (connectivity, silhouette width, Dunn 
index) then run with annotation lists to generate measures of biological
homogeneity.

See http://aging.brain-map.org/overview/home for more details about the data.

R 3.4.1
'
##

setwd('..')

'
Install packages & load libraries
'
### ~*~*~*~* RUN THIS SECTION ONLY ONCE *~*~*~*~ ###
source('https://bioconductor.org/biocLite.R')
biocLite()

# install clValid & dependencies from CRAN
install.packages(c('clValid', 'cluster', 'mclust', 'kohonen'), repo='https://CRAN.R-project.org/')

# 'mygene' - for constructing lists of functional classes
biocLite(c('Biobase', 'annotate', 'GO.db', 'mygene'))
### ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~ ###

# load libraries
library(clValid)        # Also loads dependencies
library(Biobase)        # For biological validation
library(annotate)       # For biological validation
library(GO.db)          # For biological validation
library(mygene)         # To query mygene.info for functional annotation
library(data.table)     # I/O
library(dplyr)          # Entering the TIDYVERSE!
library(gplots)         # heatmap.2

'
Load data
'
# load normalized FPKM values, sample info, DE gene lists for brain regions
fpkm_table <- readRDS(file='data/normalized_fpkm_matrix.Rds')
sample_info <- readRDS(file='data/sample_info.Rds')
brain_reg_sig_genes <- readRDS(file='data/brain_reg_sig_genes.Rds')
genes <- readRDS(file='data/genes.Rds')

# pull out individual brain region DE gene lists
hip_genes <- brain_reg_sig_genes$hip_genes
fwm_genes <- brain_reg_sig_genes$fwm_genes
pcx_genes <- brain_reg_sig_genes$pcx_genes
tcx_genes <- brain_reg_sig_genes$tcx_genes

sample_dementia_status <-  sample_info[, c('rnaseq_profile_id', 'act_demented')]

# make FPKM numeric matrix
fpkm_mat <- as.matrix(fpkm_table)
class(fpkm_mat) <- 'numeric'

# standardize
fpkm_standard_mat <- t(scale(t(fpkm_mat)))

# subset FPKM matrix by brain region
hip_samples <- sample_info$rnaseq_profile_id[which(sample_info$structure_acronym == 'HIP')]
hip_data <- fpkm_standard_mat[hip_genes, colnames(fpkm_standard_mat) %in% hip_samples]

fwm_samples <- sample_info$rnaseq_profile_id[which(sample_info$structure_acronym == 'FWM')]
fwm_data <- fpkm_standard_mat[fwm_genes, colnames(fpkm_standard_mat) %in% fwm_samples]

pcx_samples <- sample_info$rnaseq_profile_id[which(sample_info$structure_acronym == 'PCx')]
pcx_data <- fpkm_standard_mat[pcx_genes, colnames(fpkm_standard_mat) %in% pcx_samples]

tcx_samples <- sample_info$rnaseq_profile_id[which(sample_info$structure_acronym == 'TCx')]
tcx_data <- fpkm_standard_mat[tcx_genes, colnames(fpkm_standard_mat) %in% tcx_samples]

'
Gene ontology lists from mygene.info database
'
# query mygene.info ElasticSearch database to grab gene ontology info for Entrez IDs
hip_query <- queryMany(rownames(hip_data), scopes='entrezgene', fields=c('go', 'symbol'), species='human')
fwm_query <- queryMany(rownames(fwm_data), scopes='entrezgene', fields=c('go', 'symbol'), species='human')
pcx_query <- queryMany(rownames(pcx_data), scopes='entrezgene', fields=c('go', 'symbol'), species='human')
tcx_query <- queryMany(rownames(tcx_data), scopes='entrezgene', fields=c('go', 'symbol'), species='human')

# take a look at the first few entries for HIP
head(hip_query, 10)

# pull out Entrez Gene IDs
hip_entrez <- hip_query$query
fwm_entrez <- fwm_query$query
pcx_entrez <- pcx_query$query
tcx_entrez <- tcx_query$query

# take a look at a few individual gene ontology records; 
#CC for gene #1, HIP samples (Entrez Gene #25937, WWTR1)
hip_query[1, 'go.CC'][[1]]
# MF for gene #1
hip_query[1, 'go.MF'][[1]]
# BP for gene #1
hip_query[1, 'go.BP'][[1]]

# for gene #70 (Entrez Gene #23580, CDC42EP4)
hip_query[70, 'go.CC'][[1]]
hip_query[70, 'go.MF'][[1]]
hip_query[70, 'go.BP'][[1]]

# loops to extract most common functional category for each gene for three ontologies:
#    CC = cellular component
#    MF = molecular function
#    BP = biological process
### ~*~*~*~*~*~*~*~*~*~*~*~*~*~* HIP *~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~ ###
hip_cc <- NULL
hip_mf <- NULL
hip_bp <- NULL
for(i in 1:length(hip_genes)){
    # each if-else checks each 'term' list's length and, if zero, marks as 'unknown'.
    # if not zero, returns the most common term
    if(length(hip_query[i, 'go.CC'][[1]]$term) == 0){
        a_hip_cc <- 'unknown'
    }
    else{
        a_hip_cc <- tail(names(sort(table(hip_query[i, 'go.CC'][[1]]$term))), 1)
    } 
    hip_cc <- c(hip_cc, a_hip_cc)
    
    if(length(hip_query[i, 'go.MF'][[1]]$term) == 0){
        a_hip_mf <- 'unknown'
    }
    else{
        a_hip_mf <- tail(names(sort(table(hip_query[i, 'go.MF'][[1]]$term))), 1)
    } 
    hip_mf <- c(hip_mf, a_hip_mf)
        
    if(length(hip_query[i, 'go.BP'][[1]]$term) == 0){
        a_hip_bp <- 'unknown'
    }
    else{
        a_hip_bp <- tail(names(sort(table(hip_query[i, 'go.BP'][[1]]$term))), 1)
    } 
    hip_bp <- c(hip_bp, a_hip_bp)
}

hip_fc <- data.frame(cbind(hip_cc, hip_mf, hip_bp))
rownames(hip_fc) <- hip_entrez
colnames(hip_fc) <- c('CC', 'MF', 'BP')

### ~*~*~*~*~*~*~*~*~*~*~*~*~*~* FWM *~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~ ###
fwm_cc <- NULL
fwm_mf <- NULL
fwm_bp <- NULL
for(i in 1:length(fwm_genes)){
    if(length(fwm_query[i, 'go.CC'][[1]]$term) == 0){
        a_fwm_cc <- 'unknown'
    }
    else{
        a_fwm_cc <- tail(names(sort(table(fwm_query[i, 'go.CC'][[1]]$term))), 1)
    } 
    fwm_cc <- c(fwm_cc, a_fwm_cc)
    
    if(length(fwm_query[i, 'go.MF'][[1]]$term) == 0){
        a_fwm_mf <- 'unknown'
    }
    else{
        a_fwm_mf <- tail(names(sort(table(fwm_query[i, 'go.MF'][[1]]$term))), 1)
    } 
    fwm_mf <- c(fwm_mf, a_fwm_mf)
        
    if(length(fwm_query[i, 'go.BP'][[1]]$term) == 0){
        a_fwm_bp <- 'unknown'
    }
    else{
        a_fwm_bp <- tail(names(sort(table(fwm_query[i, 'go.BP'][[1]]$term))), 1)
    } 
    fwm_bp <- c(fwm_bp, a_fwm_bp)
}

fwm_fc <- data.frame(cbind(fwm_cc, fwm_mf, fwm_bp))
rownames(fwm_fc) <- fwm_entrez
colnames(fwm_fc) <- c('CC', 'MF', 'BP')

### ~*~*~*~*~*~*~*~*~*~*~*~*~*~* PCx *~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~ ###
pcx_cc <- NULL
pcx_mf <- NULL
pcx_bp <- NULL
for(i in 1:length(pcx_genes)){
    if(length(pcx_query[i, 'go.CC'][[1]]$term) == 0){
        a_pcx_cc <- 'unknown'
    }
    else{
        a_pcx_cc <- tail(names(sort(table(pcx_query[i, 'go.CC'][[1]]$term))), 1)
    } 
    pcx_cc <- c(pcx_cc, a_pcx_cc)
    
    if(length(pcx_query[i, 'go.MF'][[1]]$term) == 0){
        a_pcx_mf <- 'unknown'
    }
    else{
        a_pcx_mf <- tail(names(sort(table(pcx_query[i, 'go.MF'][[1]]$term))), 1)
    } 
    pcx_mf <- c(pcx_mf, a_pcx_mf)
        
    if(length(pcx_query[i, 'go.BP'][[1]]$term) == 0){
        a_pcx_bp <- 'unknown'
    }
    else{
        a_pcx_bp <- tail(names(sort(table(pcx_query[i, 'go.BP'][[1]]$term))), 1)
    } 
    pcx_bp <- c(pcx_bp, a_pcx_bp)
}

pcx_fc <- data.frame(cbind(pcx_cc, pcx_mf, pcx_bp))
rownames(pcx_fc) <- pcx_entrez
colnames(pcx_fc) <- c('CC', 'MF', 'BP')

### ~*~*~*~*~*~*~*~*~*~*~*~*~*~* TCx *~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~ ###
tcx_cc <- NULL
tcx_mf <- NULL
tcx_bp <- NULL
for(i in 1:length(tcx_genes)){
    if(length(tcx_query[i, 'go.CC'][[1]]$term) == 0){
        a_tcx_cc <- 'unknown'
    }
    else{
        a_tcx_cc <- tail(names(sort(table(tcx_query[i, 'go.CC'][[1]]$term))), 1)
    } 
    tcx_cc <- c(tcx_cc, a_tcx_cc)
    
    if(length(tcx_query[i, 'go.MF'][[1]]$term) == 0){
        a_tcx_mf <- 'unknown'
    }
    else{
        a_tcx_mf <- tail(names(sort(table(tcx_query[i, 'go.MF'][[1]]$term))), 1)
    } 
    tcx_mf <- c(tcx_mf, a_tcx_mf)
        
    if(length(tcx_query[i, 'go.BP'][[1]]$term) == 0){
        a_tcx_bp <- 'unknown'
    }
    else{
        a_tcx_bp <- tail(names(sort(table(tcx_query[i, 'go.BP'][[1]]$term))), 1)
    } 
    tcx_bp <- c(tcx_bp, a_tcx_bp)
}

tcx_fc <- data.frame(cbind(tcx_cc, tcx_mf, tcx_bp))
rownames(tcx_fc) <- tcx_entrez
colnames(tcx_fc) <- c('CC', 'MF', 'BP')

# Specific gene ontology list generation
# ~*~*~*~*~*~*~* Cellular Component (CC) *~*~*~*~*~*~*~* #
hip_fc_cc <- tapply(rownames(hip_data), hip_fc$CC, c)
fwm_fc_cc <- tapply(rownames(fwm_data), fwm_fc$CC, c)
pcx_fc_cc <- tapply(rownames(pcx_data), pcx_fc$CC, c)
tcx_fc_cc <- tapply(rownames(tcx_data), tcx_fc$CC, c)

# ~*~*~*~*~*~*~* Molecular Function (MF) *~*~*~*~*~*~*~* #
hip_fc_mf <- tapply(rownames(hip_data), hip_fc$MF, c)
fwm_fc_mf <- tapply(rownames(fwm_data), fwm_fc$MF, c)
pcx_fc_mf <- tapply(rownames(pcx_data), pcx_fc$MF, c)
tcx_fc_mf <- tapply(rownames(tcx_data), tcx_fc$MF, c)

# ~*~*~*~*~*~*~* Biological Process (BP) *~*~*~*~*~*~*~* #
hip_fc_bp <- tapply(rownames(hip_data), hip_fc$BP, c)
fwm_fc_bp <- tapply(rownames(fwm_data), fwm_fc$BP, c)
pcx_fc_bp <- tapply(rownames(pcx_data), pcx_fc$BP, c)
tcx_fc_bp <- tapply(rownames(tcx_data), tcx_fc$BP, c)

'
Clustering Experiment Functions
'
# exp01: Internal metrics only experiment
exp01 <- function(data, genes, k_range, clust_algos, metric, method){
    valid01 <- clValid(data, 
                       k_range, 
                       clMethods = clust_algos, 
                       validation = 'internal',
                       maxitems = length(genes),
                       metric = metric,
                       method = method,
                       verbose = FALSE)
    return(valid01)
}

# exp02: Stability metrics added
exp02 <- function(data, genes, k_range, clust_algos, metric, method){
    valid02 <- clValid(data, 
                       k_range,
                       clMethods = clust_algos,
                       validation = c('internal', 'stability'),
                       maxitems = length(genes),
                       metric = metric,
                       method = method,
                       verbose = FALSE)
    return(valid02)
}

# exp03: Biological validation using custom FC list
exp03 <- function(data, genes, k_range, clust_algos, valid_metrics, metric, method, fc_list){
    valid03 <- clValid(data,
                       k_range,
                       clMethods = clust_algos,
                       validation = valid_metrics, 
                       maxitems = length(genes),
                       metric = metric,
                       method = method,
                       annotation = fc_list,
                       verbose = FALSE)
    return(valid03)
}

'
Clustering Experiments
'
linkages <- c('complete', 'average', 'ward')
distances <- c('correlation', 'euclidean')
algos <- c('hierarchical', 'clara', 'model')

### SERIES 1
### ~*~*~*~*~*~*~*~*~*~*~*~*~*~* HIP *~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~ ###
### Internal & biological validation metrics, all 3 ontologies
hip_results01 <- NULL
hip_conditions01 <- NULL
counter <- 1
ontologies <- c('hip_fc_cc', 'hip_fc_mf', 'hip_fc_bp')

start_hip <- Sys.time()

for(linkage in linkages){
    for(distance in distances){
        for(ontology in ontologies){
            hip_conditions01[[counter]] <- list(distance = distance,
                                                linkage = linkage,
                                                ontology = ontology)
            exp <- exp03(hip_data, hip_genes, 2:12, algos, c('internal', 'biological'), distance, linkage, eval(parse(text = ontology)))
            hip_results01[counter] <- list(optimalScores(exp))
            counter <- counter + 1
            print(paste('Finished', distance, 'distance,', linkage, 'linkage,', ontology, 'ontology experiment.'))
            flush.console()
        }
    }
}

# save results and conditions to data folder
saveRDS(hip_results01, file='data/hip_clustering_results01.Rds')
saveRDS(hip_conditions01, file='data/hip_clustering_conditions01.Rds')

stop_hip <- Sys.time()
duration_hip <- stop_hip-start_hip
print(paste('Experiment duration:', duration_hip))

### ~*~*~*~*~*~*~*~*~*~*~*~*~*~* FWM *~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~ ###
# Internal metrics only
linkages <- c('average', 'ward')
distances <- c('correlation', 'euclidean')
algos <- c('hierarchical', 'clara', 'model')

fwm_results01 <- NULL
fwm_conditions01 <- NULL
counter <- 1

start_fwm <- Sys.time()

for(linkage in linkages){
    for(distance in distances){
        fwm_conditions01[[counter]] <- list(distance = distance,
                                            linkage = linkage)
        exp <- exp01(fwm_data, fwm_genes, 8:20, algos, distance, linkage)
        fwm_results01[counter] <- list(optimalScores(exp))
        counter <- counter + 1
        print(paste('Finished', distance, 'distance,', linkage, 'linkage,', 'experiment.'))
        flush.console()
    }
}

# save results and conditions to data folder
saveRDS(fwm_results01, file='data/fwm_clustering_results01.Rds')
saveRDS(fwm_conditions01, file='data/fwm_clustering_conditions01.Rds')

stop_fwm <- Sys.time()
duration_fwm <- stop_fwm-start_fwm
print(paste('Experiment duration:', duration_fwm))

### ~*~*~*~*~*~*~*~*~*~*~*~*~*~* PCx *~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~ ###


### ~*~*~*~*~*~*~*~*~*~*~*~*~*~* TCx *~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~ ###









































##########
plot(test, measure = "Dunn", legendLoc = "topleft")

hc <- clusters(test, "hierarchical")
fc_labels <- factor(hip_fc$CC)
plot(hc, labels = fc_labels)



seven <- cutree(hc, 7)
hip_clusters <- xtabs(~hip_fc$CC + seven)

'
Heatmaps
'
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
