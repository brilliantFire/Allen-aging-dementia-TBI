##
'
Differential Gene Expression Analysis ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
30 May 2018
31 May 2018 - added corrected p-values; gets top 1000 genes; produces
              summaries of gene expression difference classification

This script identifies sets of differentially expressed genes in subsets
of the Allen Aging, Dementia, and TBI RNA-seq dataset.

See http://aging.brain-map.org/overview/home for more details about the data.

R 3.4.1
'
##

'
Load libraries
'
library(edgeR)         # DGE (also loads dependencies)

'
Load data
'
setwd('..')

# Make a function to:
#   1. Load & tidy up "*_counts_status.csv" files from ../data
#   2. Extract status row & delete from orig count matrix,
#   3. Make DGEList object,
#   4. Filters counts, &
#   5. Normalize counts.
#   6. Estimate common dispersion
#   7. Estimate trended dispersions
#   8. Estimate tagwise (AKA genewise) dispersions

# Takes the following arguments:
#   group = string; group name (hip, fwm, pcx, tcx, male, or female)
#   filter_counts, fitler_samples = the minimum number of counts in the minimum number
#                                   of samples for filtering
#   ... = any other argument from any of the functions inside
prep_count_mat <- function(group, filter_counts, filter_samples, ...){
    # 1. get file & tidy it up
    x_counts_status <- data.frame(read.csv(paste('data/', group, '_counts_status.csv', sep='')))
    rownames(x_counts_status) <- x_counts_status$X
    x_counts_status$X <- NULL
    colnames(x_counts_status) <- substring(colnames(x_counts_status), 2)
    
    # 2. Extract status row and delete from count mat
    x_status <- unlist(x_counts_status[50284,])
    x_counts <- as.matrix(x_counts_status[1:50283,])
    class(x_counts) <- 'numeric'
    x_counts <- round(x_counts)
    
    # 3. Make DGEList object
    x_dgeList <- DGEList(x_counts, group = x_status)
    
    # 4. Filter
    x_filtered <- x_dgeList[rowSums(1e+06*x_dgeList$counts/expandAsMatrix(x_dgeList$samples$lib.size, dim(x_counts))>filter_counts)>=filter_samples,]
    print(paste('Number of genes left after filtering: ', dim(x_filtered)[1]))
    
    # 5. Normalize
    x_norm <- calcNormFactors(x_filtered)
    
    # 6. Estimate common dispersion
    x_norm <- estimateGLMCommonDisp(x_norm, verbose = TRUE)
    
    # 7. Estimate trended dispersion
    x_norm <- estimateGLMTrendedDisp(x_norm)
    
    # 8. Estimate tagwise ("gene-wise") dispersions
    x_norm <- estimateGLMTagwiseDisp(x_norm)
    
    # save object to .Rda file
    saveRDS(x_norm,file=paste('data/',group,'_norm_counts_disp.Rds',sep=''))
    
    # also return it for hypothesis tests below
    return(x_norm)
}

hip_norm <- prep_count_mat('hip',2,10)
fwm_norm <- prep_count_mat('fwm',2,10)
pcx_norm <- prep_count_mat('pcx',2,10)
tcx_norm <- prep_count_mat('tcx',2,10)

male_norm <- prep_count_mat('male',2,10)
female_norm <- prep_count_mat('female',2,10)

'
Hypothesis testing/DEG hunt
'
# function to perform Fisher's exact test & save results
# Takes arguments:
#   1. norm_counts = output from prep_count_mat function above
#   2. group = hip, fwm, pcx, tcx, male, or female
#   3. disp = common, trended, tagwise, or auto (auto pick the most complex dispersion available)
get_DEGs <- function(norm_counts, group, disp, ...){
    # perform Fisher's Exact Test
    x_test <- exactTest(norm_counts, pair = c('No Dementia' , 'Dementia'), dispersion = disp)
    # extract table & write to .csv in data folder
    x_table <- topTags(x_test, n=1000)
    saveRDS(x_table, file=paste('data/',group,'_exact_test_results_',disp,'.Rds',sep=''))
    return(x_test)
}

hip_results <- get_DEGs(hip_norm, 'hip', 'auto', p.value = 0.05)
fwm_results <- get_DEGs(fwm_norm, 'fwm', 'auto', p.value = 0.05)
pcx_results <- get_DEGs(pcx_norm, 'pcx', 'auto', p.value = 0.05)
tcx_results <- get_DEGs(tcx_norm, 'tcx', 'auto', p.value = 0.05)

male_results <- get_DEGs(male_norm, 'male', 'auto', p.value = 0.05)
female_results <- get_DEGs(female_norm, 'female', 'auto', p.value = 0.05)

'
Summaries
'
# summary of classification results for expression differences
summary(decideTestsDGE(hip_results, p.value = 0.05))
summary(decideTestsDGE(fwm_results, p.value = 0.05))
summary(decideTestsDGE(pcx_results, p.value = 0.05))
summary(decideTestsDGE(tcx_results, p.value = 0.05))

summary(decideTestsDGE(male_results, p.value = 0.05))
summary(decideTestsDGE(female_results, p.value = 0.05))