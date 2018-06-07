##
'
Differential Gene Expression Analysis ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
31 May 2018

This script produces the following visualizations of differentially
expressed genes identified in the XXX script...

See http://aging.brain-map.org/overview/home for more details about the data.

R 3.4.1
'
##

'
Load libraries
'
library(edgeR)
library(gplots)

'
Load dataframes produced by DE_exact_tests.R script...
'
setwd('..')
# load dataframes with normalized counts and dispersions
load_norm <- function(group, ...){
    x_norm <- readRDS(file=paste('data/',group,'_norm_counts_disp.Rds',sep=''))
    return(x_norm)
}

load_results <- function(group, disp, ...){
    x_results <- data.frame(readRDS(file=paste('data/',group,'_exact_test_results_',disp,'.Rds',sep='')))
    return(x_results)
}

hip_norm <- load_norm('hip')
hip_results <- data.frame(load_results('hip', 'auto'))

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
#fwm_bcv_plot <- bcv_plot(fwm_norm)
#pcx_bcv_plot <- bcv_plot(pcx_norm)
#tcx_bcv_plot <- bcv_plot(tcx_norm)

#male_bcv_plot <- bcv_plot(male_norm)
#female_bcv_plot <- bcv_plot(female_norm)

'
Lists of significantly different genes & summary stats
'
# Function:
#   1. Gets the names of genes with p-values from exact test below
#      a user-inputted significance threshold
#   2. Provides:
#       a. The number of significantly DE genes
#       b. The percentage of genes that are significantly DE

get_sig_genes <- function(test_results, n_genes, alpha, ...){
    # get names of genes with p <= p_value
    x_names <- rownames(test_results)[test_results$FDR <= alpha]
    # number of significant genes
    num_genes <- length(x_names)
    # percentage of total genes
    per_total <- (length(x_names)/50283)*100
    # print results
    print(paste('Number of genes w/ FDR <= ', alpha,': ', num_genes, sep=''))
    print(paste('Percentage of total # genes: ', round(per_total, 2), sep=''))
    # return list
    return(list(x_names, num_genes, per_total))
}

hip_genes <- get_sig_genes(hip_results, 500, 0.05)

'
Count histograms for the top however-many DE genes
'
hist(hip_results[unlist(hip_genes[1:200]),'logCPM'], 
     breaks=25, 
     xlab='Log Concentration (CPM)', 
     col='darkgreen', 
     freq=FALSE, 
     main='Hippocampus: Top 200')

'
MA plots
'
plotSmear(hip_norm, 
          de.tags=hip_genes, 
          main='Hippocampus', 
          pair = c('No Dementia' , 'Dementia'),
          cex = .35,
          xlab='Log Concentration (CPM)', 
          ylab='Log Fold-Change')
abline(h = c(-1, 1), col = "darkmagenta", lwd = 3)

'
TPM heatmaps
'