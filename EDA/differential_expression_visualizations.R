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
19 Jul 2018 - Changed image output format to .tiff; added HIP vs other brain regions
              and non-hippocampus regions only Venn diagrams
21 Jul 2018 - Added xtable for pretty summary tables

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

R 3.5.1
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
library(xtable)         # pretty LaTeX tables

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
# 'data/venn_brain_region_0.02_10000.tif'
# use venn.diagram from the VennDiagram library
brain_reg_venn <- venn.diagram(list(hip_genes=hip_genes,
                                    fwm_genes=fwm_genes,
                                    pcx_genes=pcx_genes,
                                    tcx_genes=tcx_genes) , 
                               filename = NULL,
                               imagetype='tiff',
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
# 'data/venn_donor_sex_0.02_10000.tif'
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
                               main='')
grid.draw(donor_sex_venn)

# get intersection genes
venn_sex_contents <- venn(list(male_genes=male_genes,
                               female_genes=female_genes), show.plot=FALSE)
venn_sex_intersections <- attr(venn_sex_contents, 'intersections')
intersection_genes_sex <- venn_sex_intersections$`male_genes:female_genes`

# compare brain region intersection genes & donor sex intersection genes
# 'data/venn_brain_region_sex_intersections_0.02_10000.tif'
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
                                      main='')
grid.draw(donor_sex_region_venn)

# Venn for HIP versus all other brain regions
# 'data/venn_brain_regions_hip_not_hip_0.02_10000.tif'
not_hip_genes <- unique(c(fwm_genes, pcx_genes, tcx_genes))
hip_not_hip_venn <- venn.diagram(list(hip_genes=hip_genes,
                                      not_hip_genes=not_hip_genes), 
                                      filename = NULL,
                                      imagetype='tiff',
                                      fill=c('red', 'blue'),
                                      alpha=c(0.3,0.3), # transparency
                                      cex = 1, 
                                      fontfamily='sans',
                                      cat.fontfamily='sans',
                                      category.names=c('hippocampus', 'other regions'), 
                                      cat.pos=c(-52,90),
                                      cat.dist=c(0.025,-0.04),
                                      cat.cex=c(0.9,0.9),
                                      main.fontfamily='sans',
                                      main='')
grid.draw(hip_not_hip_venn)

# Venn for just not HIP regions
# 'data/venn_brain_regions_not_hip_0.02_10000.tif'
not_hip_venn <- venn.diagram(list(fwm_genes=fwm_genes,
                                  pcx_genes=pcx_genes,
                                  tcx_genes=tcx_genes), 
                                      filename = NULL,
                                      imagetype='tiff',
                                      fill=c('red', 'blue', 'green'),
                                      alpha=c(0.3,0.3,0.3), # transparency
                                      cex = 1, 
                                      fontfamily='sans',
                                      cat.fontfamily='sans',
                                      category.names=c('forebrain WM', 'parietal ctx', 'temporal ctx'), 
                                      main.fontfamily='sans',
                                      main='Dementia Intersection DEGs ~ Non-hippocampus Regions')
grid.draw(not_hip_venn)

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
              hclustfun = function(x) hclust(x,method = 'complete'),
              distfun = function(x) as.dist(1-cor(t(x))), # 1 minus the Pearson corr distance
              col=bluered(20),
              na.rm=TRUE, 
              trace = 'none', 
              dendrogram = 'both', 
              key = TRUE, 
              key.title = '',
              key.xlab = '',
              key.ylab = '',
              labRow=FALSE,
              labCol=FALSE,
              denscol = 'black',
              cexRow = 0.8,
              cexCol = 0.8,
              margins = c(1,1),
              xlab='',
              ylab='')
}

# 600x600 pixel images w/ key
tiff(filename='data/hip_heatmap_wo_key.tif', height = 600, width = 600)
hip_heatmap <- DEG_heatmap(hip_samples, intersection_genes, 237, 25)
dev.off()

tiff(filename='data/fwm_heatmap_wo_key.tif', height = 600, width = 600)
DEG_heatmap(fwm_samples, intersection_genes, 555, 25)
dev.off()

tiff(filename='data/pcx_heatmap_wo_key.tif', height = 600, width = 600)
DEG_heatmap(pcx_samples, intersection_genes, 720, 25)
dev.off()

tiff(filename='data/tcx_heatmap_wo_key.tif', height = 600, width = 600)
DEG_heatmap(tcx_samples, intersection_genes, 420, 25)
dev.off()

### 700x700 hippocampus image with key, bigger axes labels, annotation

# subset FPKM values
all_hip <- as.matrix(merged[hip_genes, colnames(merged) %in% hip_samples])
class(all_hip) <- 'numeric'
all_hip <- t(scale(t(all_hip)))

hip_genes_all_samples <- as.matrix(merged[hip_genes, ])
class(hip_genes_all_samples) <- 'numeric'
hip_genes_all_samples <- t(scale(t(hip_genes_all_samples)))


# get HIP sample dementia status
hip_status <- sample_dementia_status[sample_dementia_status$rnaseq_profile_id %in% hip_samples, ]
hip_dementia_status <- hip_status$act_demented

status <- sample_dementia_status$act_demented
status_colors_1 <- replace(status, status == 'No Dementia', 'darkmagenta')
status_colors <- replace(status_colors_1, status_colors_1 == 'Dementia', 'darkgreen')

# replace with colors for heatmap annotation
hip_colors_1 <- replace(hip_dementia_status, hip_dementia_status == 'No Dementia', 'darkmagenta')
hip_colors <- replace(hip_colors_1, hip_dementia_status == 'Dementia', 'darkgreen')

# make a tiff
tiff(filename='data/hip_heatmap_ALL_GENES.tif', height = 700, width = 700)
heatmap.2(all_hip,
          hclustfun = function(x) hclust(x,method = 'complete'),
          distfun = function(x) as.dist(1-cor(t(x))), # 1 minus the Pearson corr distance
          col=bluered(20),
          #ColSideColors = hip_colors,
          denscol='black',
          na.rm=TRUE,
          scale='none',
          trace = 'none', 
          dendrogram = 'both', 
          key = TRUE, 
          key.title = '',
          key.xlab = '',
          key.ylab = '',
          cexRow = 2,
          cexCol = 2,
          margins = c(1.5,1.5),
          xlab='',
          ylab='',
          labRow=FALSE,
          labCol=FALSE)
dev.off()

### 700x700 FOREBRAIN WHITE MATTER image with key, bigger axes labels, annotation

# subset FPKM values
all_fwm <- as.matrix(merged[fwm_genes, colnames(merged) %in% fwm_samples])
class(all_fwm) <- 'numeric'
all_fwm <- t(scale(t(all_fwm)))

# get fwm sample dementia status
fwm_status <- sample_dementia_status[sample_dementia_status$rnaseq_profile_id %in% fwm_samples, ]
fwm_dementia_status <- fwm_status$act_demented

# replace with colors for heatmap annotation
fwm_colors_1 <- replace(fwm_dementia_status, fwm_dementia_status == 'No Dementia', 'darkmagenta')
fwm_colors <- replace(fwm_colors_1, fwm_dementia_status == 'Dementia', 'darkgreen')

# make a tiff
tiff(filename='data/fwm_heatmap_ALL_GENES.tif', height = 700, width = 700)
heatmap.2(all_fwm,
          hclustfun = function(x) hclust(x,method = 'complete'),
          distfun = function(x) as.dist(1-cor(t(x))), # 1 minus the Pearson corr distance
          col=bluered(20),
          ColSideColors = fwm_colors,
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
dev.off()

### 700x700 PARIETAL CORTEX image with key, bigger axes labels, annotation

# subset FPKM values
all_pcx <- as.matrix(merged[pcx_genes, colnames(merged) %in% pcx_samples])
class(all_pcx) <- 'numeric'
all_pcx <- t(scale(t(all_pcx)))

# get pcx sample dementia status
pcx_status <- sample_dementia_status[sample_dementia_status$rnaseq_profile_id %in% pcx_samples, ]
pcx_dementia_status <- pcx_status$act_demented

# replace with colors for heatmap annotation
pcx_colors_1 <- replace(pcx_dementia_status, pcx_dementia_status == 'No Dementia', 'darkmagenta')
pcx_colors <- replace(pcx_colors_1, pcx_dementia_status == 'Dementia', 'darkgreen')

# make a tiff
tiff(filename='data/pcx_heatmap_ALL_GENES.tif', height = 700, width = 700)
heatmap.2(all_pcx,
          hclustfun = function(x) hclust(x,method = 'complete'),
          distfun = function(x) as.dist(1-cor(t(x))), # 1 minus the Pearson corr distance
          col=bluered(20),
          ColSideColors = pcx_colors,
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
dev.off()

### 700x700 TEMPORAL CORTEX image with key, bigger axes labels, annotation

# subset FPKM values
all_tcx <- as.matrix(merged[tcx_genes, colnames(merged) %in% tcx_samples])
class(all_tcx) <- 'numeric'
all_tcx <- t(scale(t(all_tcx)))

# get tcx sample dementia status
tcx_status <- sample_dementia_status[sample_dementia_status$rnaseq_profile_id %in% tcx_samples, ]
tcx_dementia_status <- tcx_status$act_demented

# replace with colors for heatmap annotation
tcx_colors_1 <- replace(tcx_dementia_status, tcx_dementia_status == 'No Dementia', 'darkmagenta')
tcx_colors <- replace(tcx_colors_1, tcx_dementia_status == 'Dementia', 'darkgreen')

# make a tiff
tiff(filename='data/tcx_heatmap_ALL_GENES.tif', height = 700, width = 700)
heatmap.2(all_tcx,
          hclustfun = function(x) hclust(x,method = 'complete'),
          distfun = function(x) as.dist(1-cor(t(x))), # 1 minus the Pearson corr distance
          col=bluered(20),
          ColSideColors = tcx_colors,
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
tiff(filename='data/expression_classification_intersection_genes.tif', height = 700, width = 700)
ggplot(data = melted_class_results, aes(x=variable, y=gene_name)) + 
  theme(axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12)) +
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
Pretty summary tables
'
# donor info table
cohort_stats <- data.frame(c('2008.4 +/- 3.8', '89.0 +/- 6.3', '4.6 +/- 1.5'),
                       c('2008.7 +/- 3.4', '89.0 +/- 6.2', '4.7 +/- 2.0'))
colnames(cohort_stats) <- c('TBI Group', 'Control Group')
rownames(cohort_stats) <- c('Year of Death', 'Age (years)', 'Post-mortem Interval (hours)')

table01 <- xtable(cohort_stats)
print(table01)

print.xtable(table01, type='latex', file='data/table01.tex')

# number of samples for each donor
num_samples <- data.frame(table(rowSums(table(sample_info$donor_id,sample_info$structure_acronym))))
colnames(num_samples) <- c('Samples', 'Donors')

table02 <- xtable(t(num_samples))
align(table02) <- 'r|cccc'
print(table02, include.colnames = FALSE, hline.after = c(1))

print.xtable(table02, type='latex', file='data/table02.tex')

# n Dementia & n No Dementia sample breakdown
table(hip_norm$samples$group)
table(fwm_norm$samples$group)
table(pcx_norm$samples$group)
table(tcx_norm$samples$group)

table(male_norm$samples$group)
table(female_norm$samples$group)


all_groups_n <- data.frame(c(51, 48, 47, 51, 75, 122),
                            c(43, 45, 44, 48, 80, 100))
colnames(all_groups_n) = c('No Dementia', 'Dementia')
rownames(all_groups_n) = c('HIP', 'FWM', 'PCx', 'TCx', 'Female', 'Male')

all_groups_n['Total'] <- rowSums(all_groups_n)

table03 <- xtable(all_groups_n)
align(table03) <- 'l|cc|c'
digits(table03) <- 0
print(table03, hline.after = c(0,4))

print.xtable(table03, type='latex', file='data/table03.tex')

# Filtering & dispersions results
n_genes <- c(dim(hip_norm)[1],
             dim(fwm_norm)[1],
             dim(pcx_norm)[1],
             dim(tcx_norm)[1],
             dim(female_norm)[1],
             dim(male_norm)[1])

common_disp <- c(round(estimateGLMCommonDisp(hip_norm)$common.dispersion, 4),
                 round(estimateGLMCommonDisp(fwm_norm)$common.dispersion, 4),
                 round(estimateGLMCommonDisp(pcx_norm)$common.dispersion, 4),
                 round(estimateGLMCommonDisp(tcx_norm)$common.dispersion, 4),
                 round(estimateGLMCommonDisp(female_norm)$common.dispersion, 4),
                 round(estimateGLMCommonDisp(male_norm)$common.dispersion, 4))

bcv <- round(sqrt(common_disp),4)

n_sig_genes <- c(length(hip_genes),
                 length(fwm_genes),
                 length(pcx_genes),
                 length(tcx_genes),
                 length(female_genes),
                 length(male_genes))

all_stats <- data.frame(cbind(n_genes, common_disp, bcv, n_sig_genes))
colnames(all_stats) = c('No. Genes After Filtering',
                        '$\\hat{\\phi}_{common}$', 
                        'BCV',
                        'No. Significant Genes ($p_{adj} \\leq 0.02$)')
rownames(all_stats) = c('HIP', 'FWM', 'PCx', 'TCx', 'Female', 'Male')

table04 <- xtable(all_stats)
align(table04) <- 'l|c|c|c|c'
digits(table04) <- c(0,0,4,4,0)
print(table04,
      hline.after = c(0,4), 
      sanitize.text.function = function(x) {x})

print.xtable(table04, type='latex', file='data/table04.tex')

# intersection gene info
int_genes <- data.frame(read.csv('data/intersection_gene_info.csv'))
colnames(int_genes) <- c('Entrez Gene ID', 'Symbol', 'Name', 'Regulation in Dementia', 'General Function')

table05 <- xtable(int_genes)
print(table05,
      hline.after = c(0),
      include.rownames = FALSE,
      scalebox = 0.7)

print.xtable(table05, type='latex', file='data/table05.tex')

'
Output files
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