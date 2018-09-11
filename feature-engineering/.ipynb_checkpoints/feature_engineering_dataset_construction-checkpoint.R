##
'
Final Feature Engineering & Dataset Construction ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
07 Sep 2018 - Script created
11 Sep 2018 - Final, unprepped dataset construction code written; saves table to 
              data folder

This script runs the final clustering strategies on DE gene expression levels
for each brain region, labels each gene with its cluster assignment, and
creates new variables based on the centroids of those clusters.

See http://aging.brain-map.org/overview/home for more details about the data.

R 3.4.1
'
##

setwd('..')

'
Load libraries
'
library(data.table)     # I/O
library(cluster)

'
Load data
'
# load normalized FPKM values, sample info, DE gene lists for brain regions
fpkm_table <- readRDS(file='data/normalized_fpkm_matrix.Rds')
sample_info <- readRDS(file='data/sample_info.Rds')
brain_reg_sig_genes <- readRDS(file='data/brain_reg_sig_genes.Rds')

# gene names and annotation lists
genes <- readRDS(file='data/genes.Rds')
hip_anno <- readRDS(file='data/hip_annotation_lists.Rds')
fwm_anno <- readRDS(file='data/fwm_annotation_lists.Rds')
pcx_anno <- readRDS(file='data/pcx_annotation_lists.Rds')
tcx_anno <- readRDS(file='data/tcx_annotation_lists.Rds')

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
HIP Gene Expression Level Cluster Features
'
# run CLARA on HIP FPKM values
hip_final_clara <- clara(hip_data, k=3, metric = 'euclidean',
                       samples=1000, sampsize=250, pamLike=TRUE)

# plots
plot(hip_final_clara, color=TRUE, shade=TRUE, 
     labels=4, col.txt='black',
     col.p='black', col.clus=c('seagreen3', 'steelblue3', 'violetred3', 'sienna3'),
     cex=2, main = 'CLARA Clusters for Hippocampus DE Genes', 
     sub='',
     ylim=c(-12,12),
     lwd=2,
     lty=3)

# extract medoids
hip_cluster_variables <- hip_final_clara$medoids

# change rownames
rownames(hip_cluster_variables) <- c('gene_cluster01.HIP', 'gene_cluster02.HIP', 'gene_cluster03.HIP')

# set colnames to donor_ids corresponding to rnaseq_profile_ids (colnames for hip_cluster_variables)
colnames(hip_cluster_variables) <- sample_info$donor_id[which(sample_info$rnaseq_profile_id %in% colnames(hip_cluster_variables))]

# transpose for merging
hip_cluster_variables <- t(hip_cluster_variables)

'
FWM Gene Expression Level Cluster Features
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
     ylim=c(-12,12),
     lwd=2,
     lty=3)

# extract medoids
fwm_cluster_variables <- fwm_final_clara$medoids

# change rownames
rownames(fwm_cluster_variables) <- c('gene_cluster01.FWM', 'gene_cluster02.FWM')

# set colnames to donor_ids corresponding to rnaseq_profile_ids (colnames for fwm_cluster_variables)
colnames(fwm_cluster_variables) <- sample_info$donor_id[which(sample_info$rnaseq_profile_id %in% colnames(fwm_cluster_variables))]

# transpose for merging
fwm_cluster_variables <- t(fwm_cluster_variables)

'
PCx Gene Expression Level Cluster Features
'
# run CLARA on PCx FPKM values
pcx_final_clara <- clara(pcx_data, k=3, metric='euclidean', 
                         samples=1000, sampsize=250)

# plots
plot(pcx_final_clara, color=TRUE, shade=TRUE, 
     labels=4, col.txt='black',
     col.p='black', col.clus=c('seagreen3', 'steelblue3', 'violetred3', 'sienna3'),
     cex=2, main = 'CLARA Clusters for Parietal Cortex DE Genes', 
     sub='',
     ylim=c(-12,12),
     lwd=2,
     lty=3)

# extract medoids
pcx_cluster_variables <- pcx_final_clara$medoids

# change rownames
rownames(pcx_cluster_variables) <- c('gene_cluster01.PCx', 'gene_cluster02.PCx', 'gene_cluster03.PCx')

# set colnames to donor_ids corresponding to rnaseq_profile_ids (colnames for pcx_cluster_variables)
colnames(pcx_cluster_variables) <- sample_info$donor_id[which(sample_info$rnaseq_profile_id %in% colnames(pcx_cluster_variables))]

# transpose for merging
pcx_cluster_variables <- t(pcx_cluster_variables)

'
TCx Gene Expression Level Cluster Features
'
# run CLARA on TCx FPKM values
tcx_final_clara <- clara(tcx_data, k=2, metric='euclidean', 
                         samples=1000, sampsize=250)

# plots
plot(tcx_final_clara, color=TRUE, shade=TRUE, 
     labels=4, col.txt='black',
     col.p='black', col.clus=c('seagreen3', 'steelblue3', 'violetred3', 'sienna3'),
     cex=2, main = 'CLARA Clusters for Temporal Cortex DE Genes', 
     sub='',
     ylim=c(-12,12),
     lwd=2,
     lty=3)

# extract medoids
tcx_cluster_variables <- tcx_final_clara$medoids

# change rownames
rownames(tcx_cluster_variables) <- c('gene_cluster01.TCx', 'gene_cluster02.TCx')

# set colnames to donor_ids corresponding to rnaseq_profile_ids (colnames for tcx_cluster_variables)
colnames(tcx_cluster_variables) <- sample_info$donor_id[which(sample_info$rnaseq_profile_id %in% colnames(tcx_cluster_variables))]

# transpose for merging
tcx_cluster_variables <- t(tcx_cluster_variables)

'
Add genetic cluster features to table with other neuropathological & molecular
measurements
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
data_plus_hip <- merge(reshaped, hip_cluster_variables, by='row.names', all=TRUE)

# clean up row names
rownames(data_plus_hip) <- data_plus_hip$Row.names
data_plus_hip$Row.names <- NULL

# add FWM genetic variables
data_plus_fwm <- merge(data_plus_hip, fwm_cluster_variables, by='row.names', all=TRUE)
rownames(data_plus_fwm) <- data_plus_fwm$Row.names
data_plus_fwm$Row.names <- NULL

# add PCx genetic variables
data_plus_pcx <- merge(data_plus_fwm, pcx_cluster_variables, by='row.names', all=TRUE)
rownames(data_plus_pcx) <- data_plus_pcx$Row.names
data_plus_pcx$Row.names <- NULL

# add TCx genetic variables
data_plus_tcx <- merge(data_plus_pcx, tcx_cluster_variables, by='row.names', all=TRUE)
rownames(data_plus_tcx) <- data_plus_tcx$Row.names
data_plus_tcx$Row.names <- NULL

'
Add demographic & medical history data
'
# grab donor information table from resource website
donor_info <- data.frame(fread('http://aging.brain-map.org/api/v2/data/query.csv?criteria=model::ApiTbiDonorDetail,rma::options[num_rows$eqall]'))

# set donor_id as row names
rownames(donor_info) <- donor_info$donor_id
donor_info$donor_id <- NULL

# drop unneeded variables 'name' & 'control_set'
cols_to_drop <- c('name', 'control_set')
donor_info <- donor_info[ , !(names(donor_info) %in% cols_to_drop)]

# create final dataset (hoorary!)
final_dataset <- merge(data_plus_tcx, donor_info, by='row.names', all=TRUE)
rownames(final_dataset) <- final_dataset$Row.names
final_dataset$Row.names <- NULL

# save as .Rds
saveRDS(final_dataset, file='data/final_dataset_unprepped.Rds')