##
'
Differential Gene Expression Analysis ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
31 May 2018
29 Jun 2018 - Finalized sections
06 Jul 2018 - Added intersection gene comparison between brain region and donor
              sex groups
09 Jul 2018 - Generalized for intersection gene lists with different numbers of
              genes; tested for FDR <= 0.02 genes
13 Jul 2018 - Outputs

This script produces the following visualizations of differentially
expressed genes identified in the differential_expression_dispersions_exact_tests.R 
script. Visualizations include:
    * Venn diagrams of genes differentially expressed in dementia by brain region
      and donor sex
    * Heatmaps of normalized FPKM values for intersection genes identified from
      Venn diagram
    * Tile visualization of intersection gene regulation direction (up/down) in dementia
    * Biological coefficient of variation (BCV) plots

Additionally, it outputs the following to the `data` directory:
    * normalized_fpkm_matrix.Rds - matrix of normalized FPKM values with
      colnames = rnaseq_profile_id & rownames = gene_(entrez)_id
    * sample_info.Rds - table with info about each of the tissue samples including
      donor information such as dementia status
    * genes.Rds - table with gene info (names, database IDs, chromosome location)
    * brain_reg_sig_genes.Rds - lists of differentially-expressed genes for each brain
      region

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

'
Load dataframes produced by differential_expression_dispersions_exact_tests.R script...
'
setwd('..')

# load dataframes with normalized counts and dispersions
load_norm <- function(group, ...){
    x_norm <- readRDS(file=paste('data/',group,'_norm_counts_disp_high_low.Rds',sep=''))
    return(x_norm)
}

load_results <- function(group, disp, ...){
    x_results <- data.frame(readRDS(file=paste('data/',group,'_exact_test_results_p.02_high_low_',disp,'.Rds',sep='')))
    return(x_results)
}

# loop loads data for all 6 groups
groups <- c('hip', 'fwm', 'pcx', 'tcx', 'male', 'female')
for (group in groups){
    print(paste('Loading', group, 'group.'))
    flush.console()
    # make variables names
    norm_name <- paste(group,'_norm',sep='')
    results_name <- paste(group, '_results',sep='')
    # Assign output of load_norm and load_results to variable names
    assign(norm_name, load_norm(group))
    assign(results_name, load_results(group, 'auto'))
}

'
Lists of names of significantly different genes
'
# Function gets the names of genes with p-values from exact test below
# a user-defined significance threshold; returns the list
# 
get_sig_genes <- function(test_results, alpha, ...){
    # get names of genes with p <= p_value
    x_names <- rownames(test_results)[test_results$FDR <= alpha]
    # return list
    return(x_names)
}

# p.value = 0.1 here returns all genes in the lists
hip_genes <- get_sig_genes(hip_results, 0.1)
fwm_genes <- get_sig_genes(fwm_results, 0.1)
pcx_genes <- get_sig_genes(pcx_results, 0.1)
tcx_genes <- get_sig_genes(tcx_results, 0.1)

male_genes <- get_sig_genes(male_results, 0.1)
female_genes <- get_sig_genes(female_results, 0.1)

# number of sig genes
length(hip_genes)
length(fwm_genes)
length(pcx_genes)
length(tcx_genes)

length(male_genes)
length(female_genes)

'
Venn diagrams
'
# 'data/venn_brain_region_0.02_10000.png'
# use venn.diagram from the VennDiagram library
brain_reg_venn <- venn.diagram(list(hip_genes=hip_genes,
                                    fwm_genes=fwm_genes,
                                    pcx_genes=pcx_genes,
                                    tcx_genes=tcx_genes) , 
                               filename = NULL,
                               imagetype='png',
                               fill=c('red', 'blue', 'green', 'purple'),
                               alpha=c(0.3,0.3,0.3,0.3), # transparency
                               cex = 0.9, 
                               fontfamily='sans',
                               cat.fontfamily='sans',
                               category.names=c('hippocampus', 'forebrain WM', 'parietal ctx', 'temporal ctx'), 
                               main.fontfamily='sans',
                               main='',
                               margin=0.1)

# use filename = NULL in the above call plus grid.draw to show the plot                         
grid.draw(brain_reg_venn)

# use venn() from the gplots library to inspect contents
venn_contents <- venn(list(hip_genes=hip_genes,
                           fwm_genes=fwm_genes,
                           pcx_genes=pcx_genes,
                           tcx_genes=tcx_genes), show.plot=FALSE)
# grab the intersections
venn_intersections <- attr(venn_contents, 'intersections')
# 29 genes in the intersection of all 4 brain regions
intersection_genes <- venn_intersections$`hip_genes:fwm_genes:pcx_genes:tcx_genes`

# same now for donor sex groups
# 'data/venn_donor_sex_0.02_10000.png'
donor_sex_venn <- venn.diagram(list(male_genes=male_genes,
                                    female_genes=female_genes), 
                               filename = NULL,
                               imagetype='tiff',
                               fill=c('red', 'blue'),
                               alpha=c(0.3,0.3), # transparency
                               cex = 1, 
                               fontfamily='sans',
                               cat.fontfamily='sans',
                               category.names=c('male', 'female'), 
                               main.fontfamily='sans',
                               main='Dementia DEGs by Donor Sex')
grid.draw(donor_sex_venn)

# get intersection genes
venn_sex_contents <- venn(list(male_genes=male_genes,
                               female_genes=female_genes), show.plot=FALSE)
venn_sex_intersections <- attr(venn_sex_contents, 'intersections')
intersection_genes_sex <- venn_sex_intersections$`male_genes:female_genes`

# compare brain region intersection genes & donor sex intersection genes
# 'data/venn_brain_region_sex_intersections_0.02_10000.png'
donor_sex_region_venn <- venn.diagram(list(intersection_genes_sex=intersection_genes_sex, 
                                           intersection_genes=intersection_genes) , 
                                      filename = NULL,
                                      imagetype='tiff',
                                      fill=c('red', 'blue'),
                                      alpha=c(0.3,0.3), # transparency
                                      cex = 1, 
                                      fontfamily='sans',
                                      cat.fontfamily='sans',
                                      category.names=c('donor sex', 'brain region'), 
                                      main.fontfamily='sans',
                                      main='Dementia Intersection DEGs ~ Donor Sex & Brain Region')
grid.draw(donor_sex_region_venn)

'
FPKM Heatmaps
'
# get sample names for each group
hip_samples <- colnames(hip_norm$counts)
fwm_samples <- colnames(fwm_norm$counts)
pcx_samples <- colnames(pcx_norm$counts)
tcx_samples <- colnames(tcx_norm$counts)

male_samples <- colnames(male_norm$counts)
female_samples <- colnames(female_norm$counts)

# get normalized FPKM data (.zip file) from Allen download page
temp <- tempfile()
download.file('http://aging.brain-map.org/api/v2/well_known_file_download/502999992', temp, mode='wb')

# get list of genes out of the .zip
table1 <- unz(temp, 'rows-genes.csv')
genes <- data.frame(read.csv(table1, sep=",", header=T))

# get normalized FPKM table
table2 <- unz(temp, 'fpkm_table_normalized.csv')
fpkm <- data.frame(read.csv(table2,sep=',',header=TRUE))
unlink(temp)

# merge FPKM matrix with gene list on gene_id
colnames(fpkm)[1] <- 'gene_id'
merged <- merge(fpkm, genes, by='gene_id', all=TRUE)

# change row names to gene_entrez_id
rownames(merged) <- merged$gene_entrez_id
merged$gene_entrez_id <- NULL

# removed old id column
merged$gene_id <- NULL

# remove gene list columns
merged <- merged[,!(names(merged) %in% c('chromosome', 'gene_symbol', 'gene_name'))]

# remove X from column/sample names
colnames(merged) <- substring(colnames(merged), 2)

# Get donor info
donor_files <- data.frame(fread('http://aging.brain-map.org/api/v2/data/query.csv?criteria=model::ApiTbiDonorDetail,rma::options[num_rows$eqall]'))

# Get TBI_data_files.csv
data_files <- data.frame(fread('http://aging.brain-map.org/data/tbi_data_files.csv'))

# merge on donor_id
sample_info <- merge(data_files, donor_files, by = 'donor_id', all = TRUE)

# subset to just rnaseq_profile_id and act_demented
sample_dementia_status <-  sample_info[, c('rnaseq_profile_id', 'act_demented')]

# Heatmap function
DEG_heatmap <- function(sample_list, gene_list, RNGseed, num_samples, ...){
    # subset dementia status list
    x_dementia_status <- sample_dementia_status[sample_dementia_status$rnaseq_profile_id %in% sample_list, ]
    # dplyr library "sample_n" function to choose num_samples random samples from each of the
    # dementia status groups (dementia/no dementia)
    set.seed(RNGseed)
    x_subsample <- x_dementia_status %>% group_by(act_demented) %>% sample_n(num_samples)
    # subset FPKM by gene_list and 
    x_fpkm <- merged[gene_list, colnames(merged) %in% x_subsample$rnaseq_profile_id]
    # convert to numeric matrix
    x_fpkm <- as.matrix(x_fpkm)
    class(x_fpkm) <- 'numeric'
    # standardize
    x_fpkm <- t(scale(t(x_fpkm)))
    # heatmap
    heatmap.2(x_fpkm,
              col=bluered(20),
              na.rm=TRUE, 
              trace = 'none', 
              dendrogram = 'both', 
              key = TRUE, 
              key.title = 'Key',
              cexRow = 0.7,
              cexCol = 0.7,
              margins = c(6,6),
              xlab='samples',
              ylab='genes')
}

# 600x600 pixel images w/ key
png(filename='data/hip_heatmap_w_key.png', height = 600, width = 600)
hip_heatmap <- DEG_heatmap(hip_samples, intersection_genes, 555, 25)
dev.off()

png(filename='data/fwm_heatmap_w_key.png', height = 600, width = 600)
DEG_heatmap(fwm_samples, intersection_genes, 555, 25)
dev.off()

png(filename='data/pcx_heatmap_w_key.png', height = 600, width = 600)
DEG_heatmap(pcx_samples, intersection_genes, 720, 25)
dev.off()

png(filename='data/tcx_heatmap_w_key.png', height = 600, width = 600)
DEG_heatmap(tcx_samples, intersection_genes, 420, 25)
dev.off()

# 700x700 hippocampus image with key, bigger axes labels
png(filename='data/hip_heatmap_w_key_BIG.png', height = 700, width = 700)
hip_heatmap <- DEG_heatmap(hip_samples, intersection_genes, 555, 25)
dev.off()

'
Exact test result classification comparison for brain region intersection genes
'
# re-run exact tests on normalized count matrices for intersection genes only
# *NOTE: these are the same conditions as before just on a subset of genes.
# Adjusted p-values will differ due to different number of tests.

# function runs exact tests on the normalized counts for subsets of genes/samples & classifies 
# them as upregulated or downregulated in dementia samples
compare_gene_subset <- function(x_norm, gene_list, ...){
    x_exact_tests <- exactTest(x_norm[rownames(x_norm) %in% gene_list, ],
                               pair = c('No Dementia', 'Dementia'),
                               dispersion = 'auto')
    x_exact_test_classifications <- decideTestsDGE(x_exact_tests, p.value = 0.02)
    return(x_exact_test_classifications)
}

hip_classifications <- compare_gene_subset(hip_norm, intersection_genes)
fwm_classifications <- compare_gene_subset(fwm_norm, intersection_genes)
pcx_classifications <- compare_gene_subset(pcx_norm, intersection_genes)
tcx_classifications <- compare_gene_subset(tcx_norm, intersection_genes)

# recursive merge to construct results dataframe
hip_fwm <- merge(as.data.frame(hip_classifications),
                 as.data.frame(fwm_classifications),
                 by='row.names',
                 all=TRUE)

rownames(hip_fwm) <- hip_fwm$Row.names
hip_fwm$Row.names <- NULL
colnames(hip_fwm) <- c('HIP','FWM')


hip_fwm_pcx <- merge(hip_fwm,
                     as.data.frame(pcx_classifications),
                     by='row.names',
                     all=TRUE)

rownames(hip_fwm_pcx) <- hip_fwm_pcx$Row.names
hip_fwm_pcx$Row.names <- NULL
colnames(hip_fwm_pcx) <- c('HIP','FWM','PCx')

all_regions <- merge(hip_fwm_pcx,
                     as.data.frame(tcx_classifications),
                     by='row.names',
                     all=TRUE)

rownames(all_regions) <- all_regions$Row.names
all_regions$Row.names <- NULL
colnames(all_regions) <- c('HIP','FWM','PCx','TCx')

# make row names a column called gene_entrez_id
all_regions$gene_entrez_id <- rownames(all_regions)
rownames(all_regions) <- NULL

# Merge with gene info
merged_results <- merge(genes, all_regions, by='gene_entrez_id', all=FALSE)

# subset merged_results to just gene_name + results
class_results <- merged_results[ , c('gene_name','HIP','FWM','PCx','TCx')]

# replace 1s with 'Up' and -1s with 'Down'
class_results <- replace(class_results, class_results == 1, 'Up')
class_results <- replace(class_results, class_results == -1, 'Down')

# use melt() from the reshape2 library to reshape the dataframe
melted_class_results <- melt(class_results, id.vars='gene_name')

# use ggplot2 to visualize the classification results
png(filename='data/expression_classification_intersection_genes.png', height = 700, width = 700)
ggplot(data = melted_class_results, aes(x=variable, y=gene_name)) + 
  geom_tile(aes(fill=value), color='lightgray', show.legend=TRUE) + 
  labs(x='brain region', y='gene', fill='expression') + 
  scale_fill_manual(values=c('blue', 'red'))
dev.off()

'
Gene IDs
'
# get the rows from the `genes` table corresponding to the 29 
# intersection genes
intersection_gene_IDs <- genes[genes$gene_entrez_id %in% intersection_genes,]
intersection_gene_IDs

'
BCVPlots
'
# Function creates plots of the log2 fold changes in expression versus BCV.
# BCV is the sqrt of the estimated negative binomial dispersion.
bcv_plot <- function(norm_counts, ...){
    temp_plot <- plotBCV(norm_counts,
                         xlab = 'Avg log2 Counts per Million',
                         ylab = 'Biological CV')
    return(temp_plot)
}

hip_bcv_plot <- bcv_plot(hip_norm)
fwm_bcv_plot <- bcv_plot(fwm_norm)
pcx_bcv_plot <- bcv_plot(pcx_norm)
tcx_bcv_plot <- bcv_plot(tcx_norm)

male_bcv_plot <- bcv_plot(male_norm)
female_bcv_plot <- bcv_plot(female_norm)

'
Outputs
'
# Normalized FPKM matrix w/ colnames = rnaseq_profile_id & rownames = gene_id (Entrez ID)
saveRDS(merged, file='data/normalized_fpkm_matrix.Rds')

# sample info including act_demented
saveRDS(sample_info, file='data/sample_info.Rds')

# lists of all the significant genes from each brain region
brain_reg_sig_genes <- list(hip_genes=hip_genes,
                            fwm_genes=fwm_genes,
                            pcx_genes=pcx_genes,
                            tcx_genes=tcx_genes)
saveRDS(brain_reg_sig_genes, file='data/brain_reg_sig_genes.Rds')

# gene info
saveRDS(genes, file='data/genes.Rds')