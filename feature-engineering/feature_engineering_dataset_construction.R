##
'
Final Feature Engineering & Dataset Construction ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
07 Sep 2018 - Script created
11 Sep 2018 - Final, unprepped dataset construction code written; saves table to 
              data folder
13 Sep 2018 - Added mygene queries for medoids/cluster members
14 Sep 2018 - Annotation word clouds

This script runs the final clustering strategies on DE gene expression levels
for each brain region, labels each gene with its cluster assignment, and
creates new variables based on the medoids of those clusters. A final dataset
is constructed consisting of:
    1. Gene expression cluster medoids for 4 brain regions
    2. Neuropathological & other molecular measurements by brain region
    3. Demographic and medical history data for each donor

The final table is saved to the data folder as an .Rds file

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
library(mygene)         # To query mygene.info for functional annotation
library(wordcloud)      # Word clouds of annotation terms/words

'
Load data
'
# load normalized FPKM values, sample info, DE gene lists for brain regions
fpkm_table <- readRDS(file='data/normalized_fpkm_matrix.Rds')
sample_info <- readRDS(file='data/sample_info.Rds')
brain_reg_sig_genes <- readRDS(file='data/brain_reg_sig_genes.Rds')

# gene names
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
Add demographic & medical history data, save final dataset
'
# grab donor information table from resource website
donor_info <- data.frame(fread('http://aging.brain-map.org/api/v2/data/query.csv?criteria=model::ApiTbiDonorDetail,rma::options[num_rows$eqall]'))

# set donor_id as row names
rownames(donor_info) <- donor_info$donor_id
donor_info$donor_id <- NULL

# drop unneeded variables 'name' & 'control_set'
cols_to_drop <- c('name', 'control_set')
donor_info <- donor_info[ , !(names(donor_info) %in% cols_to_drop)]

# create final dataset (hooray!)
final_dataset <- merge(data_plus_tcx, donor_info, by='row.names', all=TRUE)
rownames(final_dataset) <- final_dataset$Row.names
final_dataset$Row.names <- NULL

# save as .Rds & .csv
saveRDS(final_dataset, file='data/final_dataset_unprepped.Rds')
write.csv(final_dataset, file='data/final_dataset_unprepped_CSV.csv')

'
mygene.info Queries/Cluster Member Biological Process Wordclouds
'
### ~*~*~*~*~*~*~*~*~*~*~* HIP *~*~*~*~*~*~*~*~*~*~*~ ###
# create cluster assignment dataframe
hip_clusters <- data.frame(cluster_assignment=hip_final_clara$clustering)

# list of genes in each cluster
hip_genes01 <- rownames(hip_clusters)[which(hip_clusters$cluster_assignment == 1)]
hip_genes02 <- rownames(hip_clusters)[which(hip_clusters$cluster_assignment == 2)]
hip_genes03 <- rownames(hip_clusters)[which(hip_clusters$cluster_assignment == 3)]

# query mygene.info
hip_query01 <- queryMany(hip_genes01, scopes='entrezgene', fields=c('go', 'symbol'), species='human')
hip_query02 <- queryMany(hip_genes02, scopes='entrezgene', fields=c('go', 'symbol'), species='human')
hip_query03 <- queryMany(hip_genes03, scopes='entrezgene', fields=c('go', 'symbol'), species='human')

# biological process (BP) gene ontology terms for each cluster
hip_cluster01_bp_terms <- data.frame(bp_terms=unlist(lapply(hip_query01$go.BP, function(x) x$term)))
hip_cluster02_bp_terms <- data.frame(bp_terms=unlist(lapply(hip_query02$go.BP, function(x) x$term)))
hip_cluster03_bp_terms <- data.frame(bp_terms=unlist(lapply(hip_query03$go.BP, function(x) x$term)))

# make tables of BP terms for each cluster
hip_terms01_table <- table(hip_cluster01_bp_terms$bp_terms)
hip_terms02_table <- table(hip_cluster02_bp_terms$bp_terms)
hip_terms03_table <- table(hip_cluster03_bp_terms$bp_terms)

# HIP cluster 01 word cloud
png('data/hip_cluster01_bp_term_wordcloud.png', height=8500, width=8500, res=300)
wordcloud(names(hip_terms01_table),as.numeric(hip_terms01_table), 
          scale=c(4, 2),min.freq=2,max.words=100, random.order=T, 
          rot.per=.10, colors=c('black', 'blue', 'red'), 
          vfont=c('sans serif','bold'))
dev.off()

# HIP cluster 02 word cloud
png('data/hip_cluster02_bp_term_wordcloud.png', height=8500, width=8500, res=300)
wordcloud(names(hip_terms02_table),as.numeric(hip_terms02_table), 
          scale=c(4, 2),min.freq=1,max.words=100, random.order=T, 
          rot.per=.10, colors=c('black', 'blue', 'red'), vfont=c('sans serif','bold'))
dev.off()

# HIP cluster 03 word cloud
png('data/hip_cluster03_bp_term_wordcloud.png', height=8500, width=8500, res=300)
wordcloud(names(hip_terms03_table),as.numeric(hip_terms03_table), 
          scale=c(4,2),min.freq=2,max.words=100, random.order=T, 
          rot.per=.10, colors=c('black', 'blue', 'red'), vfont=c('sans serif','bold'))
dev.off()
                    
### ~*~*~*~*~*~*~*~*~*~*~* FWM *~*~*~*~*~*~*~*~*~*~*~ ###
# create cluster assignment dataframe
fwm_clusters <- data.frame(cluster_assignment=fwm_final_clara$clustering)

# list of genes in each cluster
fwm_genes01 <- rownames(fwm_clusters)[which(fwm_clusters$cluster_assignment == 1)]
fwm_genes02 <- rownames(fwm_clusters)[which(fwm_clusters$cluster_assignment == 2)]

# query mygene.info
fwm_query01 <- queryMany(fwm_genes01, scopes='entrezgene', fields=c('go', 'symbol'), species='human')
fwm_query02 <- queryMany(fwm_genes02, scopes='entrezgene', fields=c('go', 'symbol'), species='human')

# biological process (BP) gene ontology terms for each cluster
fwm_cluster01_bp_terms <- data.frame(bp_terms=unlist(lapply(fwm_query01$go.BP, function(x) x$term)))
fwm_cluster02_bp_terms <- data.frame(bp_terms=unlist(lapply(fwm_query02$go.BP, function(x) x$term)))

# make tables of BP terms for each cluster
fwm_terms01_table <- table(fwm_cluster01_bp_terms$bp_terms)
fwm_terms02_table <- table(fwm_cluster02_bp_terms$bp_terms)

# FWM cluster 01 word cloud
png('data/fwm_cluster01_bp_term_wordcloud.png', height=8500, width=8500, res=300)
wordcloud(names(fwm_terms01_table),as.numeric(fwm_terms01_table), 
          scale=c(4, 2),min.freq=2,max.words=100, random.order=T, 
          rot.per=.10, colors=c('black', 'blue', 'red'), vfont=c('sans serif','bold'))
dev.off()

# FWM cluster 02 word cloud
png('data/fwm_cluster02_bp_term_wordcloud.png', height=8500, width=8500, res=300)
wordcloud(names(fwm_terms02_table),as.numeric(fwm_terms02_table), 
          scale=c(4, 2),min.freq=2,max.words=100, random.order=T, 
          rot.per=.10, colors=c('black', 'blue', 'red'), vfont=c('sans serif','bold'))
dev.off()
                                                            
### ~*~*~*~*~*~*~*~*~*~*~* PCx *~*~*~*~*~*~*~*~*~*~*~ ###
# create cluster assignment dataframe
pcx_clusters <- data.frame(cluster_assignment=pcx_final_clara$clustering)

# list of genes in each cluster
pcx_genes01 <- rownames(pcx_clusters)[which(pcx_clusters$cluster_assignment == 1)]
pcx_genes02 <- rownames(pcx_clusters)[which(pcx_clusters$cluster_assignment == 2)]
pcx_genes03 <- rownames(pcx_clusters)[which(pcx_clusters$cluster_assignment == 3)]

# query mygene.info
pcx_query01 <- queryMany(pcx_genes01, scopes='entrezgene', fields=c('go', 'symbol'), species='human')
pcx_query02 <- queryMany(pcx_genes02, scopes='entrezgene', fields=c('go', 'symbol'), species='human')
pcx_query03 <- queryMany(pcx_genes03, scopes='entrezgene', fields=c('go', 'symbol'), species='human')

# biological process (BP) gene ontology terms for each cluster
pcx_cluster01_bp_terms <- data.frame(bp_terms=unlist(lapply(pcx_query01$go.BP, function(x) x$term)))
pcx_cluster02_bp_terms <- data.frame(bp_terms=unlist(lapply(pcx_query02$go.BP, function(x) x$term)))
pcx_cluster03_bp_terms <- data.frame(bp_terms=unlist(lapply(pcx_query03$go.BP, function(x) x$term)))

# make tables of BP terms for each cluster
pcx_terms01_table <- table(pcx_cluster01_bp_terms$bp_terms)
pcx_terms02_table <- table(pcx_cluster02_bp_terms$bp_terms)
pcx_terms03_table <- table(pcx_cluster03_bp_terms$bp_terms)

# PCx cluster 01 word cloud
png('data/pcx_cluster01_bp_term_wordcloud.png', height=8500, width=8500, res=300)
wordcloud(names(pcx_terms01_table),as.numeric(pcx_terms01_table), 
          scale=c(4, 2),min.freq=2,max.words=100, random.order=T, 
          rot.per=.10, colors=c('black', 'blue', 'red'), vfont=c('sans serif','bold'))
dev.off()

# PCx cluster 02 word cloud
png('data/pcx_cluster02_bp_term_wordcloud.png', height=8500, width=8500, res=300)
wordcloud(names(pcx_terms02_table),as.numeric(pcx_terms02_table), 
          scale=c(4, 2),min.freq=1,max.words=100, random.order=T, 
          rot.per=.10, colors=c('black', 'blue', 'red'), vfont=c('sans serif','bold'))
dev.off()

# PCx cluster 03 word cloud
png('data/pcx_cluster03_bp_term_wordcloud.png', height=8500, width=8500, res=300)
wordcloud(names(pcx_terms03_table),as.numeric(pcx_terms03_table), 
          scale=c(4,2),min.freq=2,max.words=100, random.order=T, 
          rot.per=.10, colors=c('black', 'blue', 'red'), vfont=c('sans serif','bold'))
dev.off()
                                                            
### ~*~*~*~*~*~*~*~*~*~*~* TCx *~*~*~*~*~*~*~*~*~*~*~ ###
# create cluster assignment dataframe
tcx_clusters <- data.frame(cluster_assignment=tcx_final_clara$clustering)

# list of genes in each cluster
tcx_genes01 <- rownames(tcx_clusters)[which(tcx_clusters$cluster_assignment == 1)]
tcx_genes02 <- rownames(tcx_clusters)[which(tcx_clusters$cluster_assignment == 2)]

# query mygene.info
tcx_query01 <- queryMany(tcx_genes01, scopes='entrezgene', fields=c('go', 'symbol'), species='human')
tcx_query02 <- queryMany(tcx_genes02, scopes='entrezgene', fields=c('go', 'symbol'), species='human')

# biological process (BP) gene ontology terms for each cluster
tcx_cluster01_bp_terms <- data.frame(bp_terms=unlist(lapply(tcx_query01$go.BP, function(x) x$term)))
tcx_cluster02_bp_terms <- data.frame(bp_terms=unlist(lapply(tcx_query02$go.BP, function(x) x$term)))

# make tables of BP terms for each cluster
tcx_terms01_table <- table(tcx_cluster01_bp_terms$bp_terms)
tcx_terms02_table <- table(tcx_cluster02_bp_terms$bp_terms)

# tcx cluster 01 word cloud
png('data/tcx_cluster01_bp_term_wordcloud.png', height=8500, width=8500, res=300)
wordcloud(names(tcx_terms01_table),as.numeric(tcx_terms01_table), 
          scale=c(4, 2),min.freq=2,max.words=100, random.order=T, 
          rot.per=.10, colors=c('black', 'blue', 'red'), vfont=c('sans serif','bold'))
dev.off()

# tcx cluster 02 word cloud
png('data/tcx_cluster02_bp_term_wordcloud.png', height=8500, width=8500, res=300)
wordcloud(names(tcx_terms02_table),as.numeric(tcx_terms02_table), 
          scale=c(4, 2),min.freq=2,max.words=100, random.order=T, 
          rot.per=.10, colors=c('black', 'blue', 'red'), vfont=c('sans serif','bold'))
dev.off()