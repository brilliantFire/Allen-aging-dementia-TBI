##
'
Final Feature Engineering & Dataset Construction ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
07 Sep 2018 - Script created

This script runs the final clustering strategies on DE gene expression levels
for each brain region, labels each gene with its cluster assignment, and
creates new variables based on the centroids of those clusters.

See http://aging.brain-map.org/overview/home for more details about the data.

R 3.4.1
'
##

setwd('..')

'
Install packages & load libraries
'
# install ggfortify package for plots
install.packages('ggfortify', repo='https://CRAN.R-project.org/')


library(data.table)     # I/O
library(cluster)
library(ggfortify)
library(ggplot2)

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
HIP Features
'
# run CLARA on HIP FPKM values
hip_final_pam <- clara(hip_data, k=3, metric = 'euclidean',
                       samples=1000, sampsize=250, pamLike=TRUE)

# plots
plot(hip_final_pam, color=TRUE, shade=TRUE, 
     labels=4, col.txt='black',
     col.p='black', col.clus=c('seagreen3', 'steelblue3', 'violetred3', 'sienna3'),
     cex=2, main = 'PAM Clusters for Hippocampus DE Genes', 
     sub='',
     ylim=c(-12,12))

# extract cluster assignments
hip_clusters <- data.frame(hip_final_pam$clustering)
colnames(hip_clusters) <- 'cluster_assignment'

# generate variables from clusterwise means
hip_cluster_variables <- NULL
for(cluster_number in 1:max(hip_clusters$cluster_assignment)){
    clusterwise_means <- colMeans(hip_data[which(hip_clusters$cluster_assignment == cluster_number),])
    hip_cluster_variables <- rbind(hip_cluster_variables, clusterwise_means)
}

# change rownames
rownames(hip_cluster_variables) <- c('gene_cluster01.HIP', 'gene_cluster02.HIP', 'gene_cluster03.HIP')

# set colnames to donor_ids corresponding to rnaseq_profile_ids (colnames for hip_cluster_variables)
colnames(hip_cluster_variables) <- sample_info$donor_id[which(sample_info$rnaseq_profile_id %in% colnames(hip_cluster_variables))]

# transpose for merging
hip_cluster_variables <- t(hip_cluster_variables)

'
FWM Features
'
# run CLARA on FWM FPKM values
fwm_final_clara <- clara(fwm_data, k=2, metric = 'euclidean',
                         samples=1000, sampsize=250, pamLike=TRUE)

# plots
plot(fwm_final_clara, color=TRUE, shade=TRUE, 
     labels=4, col.txt='black',
     col.p='black', col.clus=c('seagreen3', 'steelblue3', 'violetred3', 'sienna3'),
     cex=2, main = 'CLARA Clusters for Forebrain White Matter DE Genes', 
     sub='',
     ylim=c(-12,12))

# extract cluster assignments
fwm_clusters <- data.frame(fwm_final_pam$clustering)
colnames(fwm_clusters) <- 'cluster_assignment'

# generate variables from clusterwise means
fwm_cluster_variables <- NULL
for(cluster_number in 1:max(fwm_clusters$cluster_assignment)){
    clusterwise_means <- colMeans(fwm_data[which(fwm_clusters$cluster_assignment == cluster_number),])
    fwm_cluster_variables <- rbind(fwm_cluster_variables, clusterwise_means)
}

# change rownames
rownames(fwm_cluster_variables) <- c('gene_cluster01.FWM', 'gene_cluster02.FWM')

# set colnames to donor_ids corresponding to rnaseq_profile_ids (colnames for fwm_cluster_variables)
colnames(fwm_cluster_variables) <- sample_info$donor_id[which(sample_info$rnaseq_profile_id %in% colnames(fwm_cluster_variables))]

# transpose for merging
fwm_cluster_variables <- t(fwm_cluster_variables)

'
PCx Features
'
# run CLARA on PCx FPKM values
pcx_final_clara <- clara(pcx_data, k=3, metric='euclidean', 
                         samples=1000, sampsize=250)

# plots
plot(pcx_final_clara, color=TRUE, shade=TRUE, 
     labels=4, col.txt='black',
     col.p='black', col.clus=c('seagreen3', 'steelblue3', 'violetred3', 'sienna3'),
     cex=2, main = 'CLARA Clusters for Parietal Cortex DE Genes', 
     #sub='',
     ylim=c(-12,12))

# extract cluster assignments
pcx_clusters <- data.frame(pcx_final_pam$clustering)
colnames(pcx_clusters) <- 'cluster_assignment'

# generate variables from clusterwise means
pcx_cluster_variables <- NULL
for(cluster_number in 1:max(pcx_clusters$cluster_assignment)){
    clusterwise_means <- colMeans(pcx_data[which(pcx_clusters$cluster_assignment == cluster_number),])
    pcx_cluster_variables <- rbind(pcx_cluster_variables, clusterwise_means)
}

# change rownames
rownames(pcx_cluster_variables) <- c('gene_cluster01.PCx', 'gene_cluster02.PCx')

# set colnames to donor_ids corresponding to rnaseq_profile_ids (colnames for pcx_cluster_variables)
colnames(pcx_cluster_variables) <- sample_info$donor_id[which(sample_info$rnaseq_profile_id %in% colnames(pcx_cluster_variables))]

# transpose for merging
pcx_cluster_variables <- t(pcx_cluster_variables)

'
Construct full dataset with new genetic features
'
# Load Luminex protein, immunohistochemistry, and isoprostane quants
neuropath_data <- data.frame(fread('http://aging.brain-map.org/api/v2/data/query.csv?criteria=model::ApiTbiDonorMetric,rma::options[num_rows$eqall]'))

# drop donor_name and structure_id (not needed)
cols_to_remove <- c('donor_name', 'structure_id')
neuropath_data <- neuropath_data[ , !(names(neuropath_data) %in% cols_to_remove)]

# reshape to wide format, donor_id down the side, variables by brain region across the top
reshaped <- reshape(neuropath_data, direction = 'wide', 
                    idvar = 'donor_id', timevar = 'structure_acronym')

# set donor_id as row names
rownames(reshaped) <- reshaped$donor_id
reshaped$donor_id <- NULL

# add HIP genetic variables
data_plus_hip <- merge(reshaped, hip_cluster_variables, by='row.names')

# clean up row names
rownames(data_plus_hip) <- data_plus_hip$Row.names
data_plus_hip$Row.names <- NULL

# all FWM genetic variables
data_plus_fwm <- merge(data_plus_hip, fwm_cluster_variables, by='row.names')
rownames(data_plus_fwm) <- data_plus_fwm$Row.names
data_plus_fwm$Row.names <- NULL

