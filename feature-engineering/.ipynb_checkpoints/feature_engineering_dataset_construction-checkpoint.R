##
'
Final Feature Engineering & Dataset Construction ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
07 Sep 2018 - Script created
11 Sep 2018 - Final, unprepped dataset construction code written; saves table to 
              data folder
13 Sep 2018 - Added mygene queries for medoids/cluster members
14 Sep 2018 - Annotation word clouds
18 Sep 2018 - Added pretty bivariate & silhouette plots for thesis figures
30 Sep 2018 - Added xtables for pretty LaTeX tables
07 Oct 2018 - Added observation weights to final_dataset

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
library(xtable)         # pretty LaTeX tables

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
     #sub='',
     ylim=c(-12,12),
     lwd=2,
     lty=3)

# pretty clusplot
png('data/hip_pretty_bivariate_plot.png', height=500, width=500)
par(mar=c(rep(4.5,4)))
clusplot(hip_final_clara, color=TRUE, shade=TRUE, 
     labels=4, col.txt='black',
     col.p='black', col.clus=c('seagreen3', 'steelblue3', 'violetred3', 'sienna3'),
     cex=2, cex.axis=1.5, cex.lab=1.5, main = '', 
     sub='',
     ylim=c(-12,12),
     lwd=2,
     lty=3)
dev.off()

# pretty silhouette plot
png('data/hip_pretty_silhouette_plot.png', height=500, width=500)
plot(silhouette(hip_final_clara, full=FALSE),
           col=c('seagreen3', 'steelblue3', 'violetred3'),
           cex=3, cex.axis=1.5, cex.lab=1.5, main = '', 
           sub='',
           xlab='',
           lwd=2,
           lty=1, do.col.sort=TRUE, do.clus.stat=FALSE, do.n.k=FALSE)
grid(col='darkgray', lty=2)
dev.off()

# extract medoids
hip_cluster_variables <- hip_final_clara$medoids

# query mygene.info for medoids
hip_medoid_query <- queryMany(rownames(hip_cluster_variables), scopes='entrezgene', fields=c('go', 'symbol'), species='human')

# get medoid GOs
hip_medoid_bp_terms <- lapply(hip_medoid_query$go.BP, function(x) x$term)
hip_medoid_cc_terms <- lapply(hip_medoid_query$go.CC, function(x) x$term)
hip_medoid_mf_terms <- lapply(hip_medoid_query$go.MF, function(x) x$term)                             

# get medoid gene info
hip_medoid_gene_info <- genes[genes$gene_entrez_id %in% rownames(hip_cluster_variables), ]

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
     #sub='',
     ylim=c(-12,12),
     lwd=2,
     lty=3)

# pretty clusplot
png('data/fwm_pretty_bivariate_plot.png', height=500, width=500)
par(mar=c(rep(4.5,4)))
clusplot(fwm_final_clara, color=TRUE, shade=TRUE, 
     labels=4, col.txt='black',
     col.p='black', col.clus=c('seagreen3', 'steelblue3', 'violetred3', 'sienna3'),
     cex=2, cex.axis=1.5, cex.lab=1.5, main = '', 
     sub='',
     ylim=c(-12,12),
     lwd=2,
     lty=3)
dev.off()

# pretty silhouette plot
png('data/fwm_pretty_silhouette_plot.png', height=500, width=500)
plot(silhouette(fwm_final_clara, full=FALSE),
           col=c('steelblue3', 'seagreen3'),
           cex=3, cex.axis=1.5, cex.lab=1.5, main = '', 
           sub='',
           xlab='',
           lwd=2,
           lty=1, do.col.sort=TRUE, do.clus.stat=FALSE, do.n.k=FALSE)
grid(col='darkgray', lty=2)
dev.off()
                              
# extract medoids
fwm_cluster_variables <- fwm_final_clara$medoids

# query mygene.info for medoids
fwm_medoid_query <- queryMany(rownames(fwm_cluster_variables), scopes='entrezgene', fields=c('go', 'symbol'), species='human')

# get medoid GOs
fwm_medoid_bp_terms <- lapply(fwm_medoid_query$go.BP, function(x) x$term)
fwm_medoid_cc_terms <- lapply(fwm_medoid_query$go.CC, function(x) x$term)
fwm_medoid_mf_terms <- lapply(fwm_medoid_query$go.MF, function(x) x$term)                             

# get medoid gene info
fwm_medoid_gene_info <- genes[genes$gene_entrez_id %in% rownames(fwm_cluster_variables), ]

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
     #sub='',
     ylim=c(-12,12),
     lwd=2,
     lty=3)

# pretty clusplot
png('data/pcx_pretty_bivariate_plot.png', height=500, width=500)
par(mar=c(rep(4.5,4)))
clusplot(pcx_final_clara, color=TRUE, shade=TRUE, 
     labels=4, col.txt='black',
     col.p='black', col.clus=c('seagreen3', 'steelblue3', 'violetred3'),
     cex=2, cex.axis=1.5, cex.lab=1.5, main = '', 
     sub='',
     ylim=c(-12,12),
     lwd=2,
     lty=3)
dev.off()

# pretty silhouette plot
png('data/pcx_pretty_silhouette_plot.png', height=500, width=500)
plot(silhouette(pcx_final_clara, full=FALSE),
           col=c('steelblue3', 'violetred3', 'seagreen3'),
           cex=3, cex.axis=1.5, cex.lab=1.5, main = '', 
           sub='',
           xlab='',
           lwd=2,
           lty=1, do.col.sort=TRUE, do.clus.stat=FALSE, do.n.k=FALSE)
grid(col='darkgray', lty=2)
dev.off()
                              
# extract medoids
pcx_cluster_variables <- pcx_final_clara$medoids

# query mygene.info for medoids
pcx_medoid_query <- queryMany(rownames(pcx_cluster_variables), scopes='entrezgene', fields=c('go', 'symbol'), species='human')

# get medoid GOs
pcx_medoid_bp_terms <- lapply(pcx_medoid_query$go.BP, function(x) x$term)
pcx_medoid_cc_terms <- lapply(pcx_medoid_query$go.CC, function(x) x$term)
pcx_medoid_mf_terms <- lapply(pcx_medoid_query$go.MF, function(x) x$term)                             

# get medoid gene info
pcx_medoid_gene_info <- genes[genes$gene_entrez_id %in% rownames(pcx_cluster_variables), ]

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
     #sub='',
     ylim=c(-12,12),
     lwd=2,
     lty=3)

# pretty clusplot
png('data/tcx_pretty_bivariate_plot.png', height=500, width=500)
par(mar=c(rep(4.5,4)))
clusplot(tcx_final_clara, color=TRUE, shade=TRUE, 
     labels=4, col.txt='black',
     col.p='black', col.clus=c('seagreen3', 'steelblue3', 'violetred3', 'sienna3'),
     cex=2, cex.axis=1.5, cex.lab=1.5, main = '', 
     sub='',
     ylim=c(-12,12),
     lwd=2,
     lty=3)
dev.off()                              

# pretty silhouette plot
png('data/tcx_pretty_silhouette_plot.png', height=500, width=500)
plot(silhouette(tcx_final_clara, full=FALSE),
           col=c('seagreen3', 'steelblue3'),
           cex=3, cex.axis=1.5, cex.lab=1.5, main = '', 
           sub='',
           xlab='',
           lwd=2,
           lty=1, do.col.sort=TRUE, do.clus.stat=FALSE, do.n.k=FALSE)
grid(col='darkgray', lty=2)
dev.off()
                              
# extract medoids
tcx_cluster_variables <- tcx_final_clara$medoids

# query mygene.info for medoids
tcx_medoid_query <- queryMany(rownames(tcx_cluster_variables), scopes='entrezgene', fields=c('go', 'symbol'), species='human')

# get medoid GOs
tcx_medoid_bp_terms <- lapply(tcx_medoid_query$go.BP, function(x) x$term)
tcx_medoid_cc_terms <- lapply(tcx_medoid_query$go.CC, function(x) x$term)
tcx_medoid_mf_terms <- lapply(tcx_medoid_query$go.MF, function(x) x$term)                             

# get medoid gene info
tcx_medoid_gene_info <- genes[genes$gene_entrez_id %in% rownames(tcx_cluster_variables), ]
                              
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
                              
#### Get local copy
neuropath_data <- read.csv('data/ProteinAndPathologyQuantifications.csv')
                              
#### Get local copy
obs_weights <- read.csv('data/group_weights.csv')
colnames(obs_weights) <- c('donor_name', 'obs_weight')

# merge neuropath_data & obs_weights on donor_name
neuropath_data <- merge(neuropath_data, obs_weights, by='donor_name', all.x=TRUE)
                              
# drop structure_id (not needed)
neuropath_data$structure_id <- NULL

# reshape to wide format, donor_id down the side, variables by brain region across the top
reshaped <- reshape(neuropath_data, direction = 'wide', 
                    idvar = 'donor_id', timevar = 'structure_acronym')
                              
# add mean of obs_weight across structures
reshaped$obs_weight <- rowMeans(cbind(reshaped$obs_weight.HIP,
                                      reshaped$obs_weight.FWM,
                                      reshaped$obs_weight.PCx,
                                      reshaped$obs_weight.TCx), na.rm=TRUE)
                              
# remove all but one donor_name variable & change name
cols_to_remove <- c('donor_name.FWM', 'donor_name.PCx', 'donor_name.TCx', 'donor_name.HIP',
                    'obs_weight.FWM', 'obs_weight.PCx', 'obs_weight.TCx', 'obs_weight.HIP')
reshaped <- reshaped[ , !(names(reshaped) %in% cols_to_remove)]
                              
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

#### Get local copy
donor_info <- read.csv('data/DonorInformation.csv')
                              
# set donor_id as row names
rownames(donor_info) <- donor_info$donor_id
donor_info$donor_id <- NULL
                              
# drop unneeded variables 'control_set' and 'name'
donor_info$control_set <- NULL
donor_info$name <- NULL

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
                              
# medoid cluster assignments
hip_medoid_assignments <- hip_final_clara$clustering[rownames(hip_clusters) %in% rownames(hip_final_clara$medoids)]

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
                                                            
# medoid cluster assignments
fwm_medoid_assignments <- fwm_final_clara$clustering[rownames(fwm_clusters) %in% rownames(fwm_final_clara$medoids)]

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
                                                            
# medoid cluster assignments
pcx_medoid_assignments <- pcx_final_clara$clustering[rownames(pcx_clusters) %in% rownames(pcx_final_clara$medoids)]

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
                                                            
# medoid cluster assignments
tcx_medoid_assignments <- tcx_final_clara$clustering[rownames(tcx_clusters) %in% rownames(tcx_final_clara$medoids)]

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
                                                            
'
Pretty LaTeX tables for thesis
'
# Originating clusters
orig_clusters <- c('HIP 1', 'HIP 2', 'HIP 3',
                   'FWM 1', 'FWM 2',
                   'PCx 1', 'PCx 2', 'PCx 3',
                   'TCx 1', 'TCx 2')
                                                            
# Entrez IDs
entrez_IDs <- c(rownames(hip_final_clara$medoids),
                rownames(fwm_final_clara$medoids),
                rownames(pcx_final_clara$medoids),
                rownames(tcx_final_clara$medoids))
                                                            
# get gene info
gene_info <- genes[match(entrez_IDs, genes$gene_entrez_id),]

# new variable names
new_vars <- c('gene_cluster01.HIP', 'gene_cluster02.HIP', 'gene_cluster03.HIP',
              'gene_cluster01.FWM', 'gene_cluster02.FWM',
              'gene_cluster01.PCx', 'gene_cluster02.PCx', 'gene_cluster03.PCx',
              'gene_cluster01.TCx', 'gene_cluster02.TCx')

# combine & clean u
medoid_df <- data.frame(cbind(orig_clusters, gene_info, new_vars))
medoid_df <- medoid_df[ , !names(medoid_df) %in% c('gene_id')]
                                                            
medoid_df <- medoid_df[ , c('orig_clusters',
                            'gene_entrez_id',
                            'gene_symbol',
                            'gene_name',
                            'chromosome',
                            'new_vars')]
rownames(medoid_df) <- NULL
colnames(medoid_df) <- c('Originating Cluster',
                         'Entrez Gene ID',
                         'Symbol',
                         'Name',
                         'Chromosome Location',
                         'New Variable Name')
                                                            
table06 <- xtable(medoid_df)
align(table06) <- 'l|c|c|c|l|c|l'
digits(table06) <- 0
print(table06, hline.after = c(-1,0,3,5,8,10), include.rownames=FALSE)

print.xtable(table06, type='latex', file='data/table06.tex')

                                                            
                                                            
                                                            
                                                            
                                                            
                                                            
                                                            
                                                            
                                                            
                                                            