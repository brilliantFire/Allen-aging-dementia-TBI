##
'
Differential Gene Expression Analysis ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
30 May 2018

This script produces the following visualizations of differentially
expressed genes identified in the XXX script...

See http://aging.brain-map.org/overview/home for more details about the data.

R 3.4.1
'
##
'
Load dataframes produced by DE_exact_tests.R script
'
setwd('..')
# load dataframes with normalized counts and dispersions
load_norm <- function(group, ...){
    x_norm <- readRDS(file=paste('data/',group,'_norm_counts_disp.Rds',sep=''))
    return(x_norm)
}

load_results <- function(group, disp, ...){
    x_results <- readRDS(file=paste('data/',group,'_exact_test_results_',disp,'.Rds',sep=''))
    return(x_results)
}

hip_norm <- load_norm('hip')
hip_results <- load_results('hip', 'auto')

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
Summary statistics for top X DE genes
'
# Function:
#   1. Gets the names of genes with p-values from exact test below
#      a user-inputted significance threshold
#   2. Provides:
#       a. The number of significantly DE genes
#       b. The percentage of genes that are significantly DE
#       c. A summary of the number of upregulated/downregulated/NS
#          ('not significant') genes

# sort by PValue
hip_common_table <- hip_common_table[order(hip_common_table$PValue),]
hip_tagwise_table <- hip_tagwise_table[order(hip_tagwise_table$PValue),]

# Names of genes where the p-value is significant at the 0.01 level
names_hip_common <- unlist(rownames(hip_common_table)[hip_common_table$PValue <= 0.01])
names_hip_tagwise <- rownames(hip_tagwise_table)[hip_tagwise_table$PValue <= 0.01]

# number of DEGs
length(names_hip_common)
length(names_hip_tagwise)

# decideTestsDGE classifies the fold change metrics in the table as up, down, or NS;
# uses correction for repeated measures that controls the false discovery rate 
# (expected proportion of false discoveries among rejected hypotheses).
hip_common_decide <- decideTestsDGE(hip_common_test, p.value = 0.01)


'
Count histograms for the top 100 DE genes
'
hist(results_common[diff_genes_common[1:100],"logCPM"], 
     breaks=25, 
     xlab="Log Concentration (CPM)", 
     col="darkgreen", 
     freq=FALSE, 
     main="Common: Top 100")