##
'
Differential Gene Expression Analysis ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
05 May 2018

Performs differential gene expression analysis using the DESeq2 package...

See http://aging.brain-map.org/overview/home for more details about the data.

R 3.4.3
'
##
'
Get DGE packages from bioconductor.org & load libraries
'
# get DGE packages from bioconductor this gets bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite()

# Now get packages
biocLite(c('edgeR', 'DESeq2', 'limma'))

# load libraries
library(data.table)    # i/o
library(edgeR)         # DGE
library(gplots)        # heatmap.2

'
Prep gene-sample matrix
'
# Add data load from csv as matrix
# raw_read_counts <- data.frame(read.csv('C:/Users/rlvis/Documents/NU_MSPA/thesis/data/raw_read_counts.csv')) # desktop
raw_read_counts <- data.frame(read.csv('C:/Users/Rebecca/Documents/NU_MSPA/thesis/data/raw_read_counts.csv')) # laptop

# make 'gene_id' row names
rownames(raw_read_counts) <- raw_read_counts$gene_id
raw_read_counts$gene_id <- NULL
head(raw_read_counts,10)

# Remove X from column names
colnames(raw_read_counts) <- substring(colnames(raw_read_counts), 2)
head(raw_read_counts,10)

'
Data overview
'
# Size of the matrix
dim(raw_read_counts)   # 50283x378

# sum down columns to get library sizes (in millions of reads)
lib_sizes = (colSums(raw_read_counts[,-1])) / 1e06
summary(lib_sizes)
barplot(lib_sizes)

# genes with low counts
table(rowSums(raw_read_counts[,-1]))[1:40]

'
Get donor & brain sample info; merge & match to gene matrix column names/rnaseq_profile_id; get dementia status
'
# Get donor info
donor_files <- data.frame(fread('http://aging.brain-map.org/api/v2/data/query.csv?criteria=model::ApiTbiDonorDetail,rma::options[num_rows$eqall]'))

# Get TBI_data_files.csv
data_files <- data.frame(fread('http://aging.brain-map.org/data/tbi_data_files.csv'))

# merge on donor_id
sample_info <- merge(data_files, donor_files, by = 'donor_id', all = TRUE)
head(sample_info,10)

# add gene matrix colnames as a column to sample_info; subtract and sum to check order
sample_info['gene_matrix_cols'] = colnames(raw_read_counts)
sample_info$order_check <- (as.numeric(sample_info$gene_matrix_cols) - as.numeric(sample_info$rnaseq_profile_id))
sum(sample_info$order_check) # 0 - we can proceed

# get target variable act_demented
dementia_status <- t(sample_info['act_demented'])

'
Constrct edgeR object & filter out low count reads
'
# puts counts into 'counts' & group into 'samples'
DGE_list <- DGEList(raw_read_counts, group = dementia_status)

# number of genes with zero counts in every sample
DGE_list$all.zeros

# filter out any genes with fewer than 1 read/million in at least 5 samples
filtered_list <- DGE_list[rowSums(1e+06 * DGE_list$counts/expandAsMatrix(DGE_list$samples$lib.size, dim(raw_read_counts)) > 1) >= 5, ]
dim(filtered_list)

'
Calculate normalization factors & effective library sizes
'
# calculate normalization factors
filtered_list <- calcNormFactors(filtered_list)
filtered_list$samples

# Now multiply by actual library size to get the effective library size
effective_lib_sizes <- filtered_list$samples$lib.size * filtered_list$samples$norm.factors
summary(effective_lib_sizes)

'
Multi-dimensional scaling plot
'
# Calculates a similarity measure & projects into 2D
plotMDS(filtered_list, method = "bcv", main = "MDS Plot for Count Data", labels = as.numeric(colnames(filtered_list$counts)))

'
Estimate the common & tagwise dispersions
'
# estimate the common dispersion (negative binomial model)
filtered_list <- estimateCommonDisp(filtered_list)
names(filtered_list)

# The estimate is in common.dispersion
filtered_list$common.dispersion

# use the common dispersion to calculate individual 'tagwise' dispersions for each gene
# prior df
filtered_list <- estimateTagwiseDisp(filtered_list, prior.df = 40)
names(filtered_list)
summary(filtered_list$tagwise.dispersion)

'
Mean-variance Biological coefficient of variation (BCV) plots
'
# Plots 4 different variance estimates: 
    # 1. The raw variances of the counts (grey dots);
    # 2. The variances using the tagwise dispersions (light blue dots);
    # 3. The variances using the common dispersion (solid blue line); and, 
    # 4. The variance = mean a.k.a. poisson variance (solid black line).

mean_var_plot <- plotMeanVar(filtered_list, 
                           show.raw.vars=TRUE,
                           show.tagwise.vars=TRUE,
                           show.binned.common.disp.vars=FALSE,
                           show.ave.raw.vars=FALSE,
                           NBline = TRUE,
                           nbins = 100,
                           pch = 16,
                           xlab ="Mean Expression (log10)",
                           ylab = "Expression Variance (log10)",
                           main = "Mean-Variance Plot" )

plotBCV(filtered_list)

'
Perform pairwise hypothesis testing for differential expression between groups 
'
# Fisher's exact test is used to test whether the expression levels for a gene 
# are different between the two groups (dementia & no dementia).

# Try three different variance estimates in the test:
    # 1. The common dispersion (same estimated value for each gene)
    # 2. The tagwise, or gene-wise, dispersions (different values for each gene)
    # 3. The Poission distribution dispersion which is = 0
# Also use the dispersion = "auto" option, which uses the most complex dispersions

diff_ex_common <- exactTest(filtered_list, pair = c("No Dementia" , "Dementia"), dispersion = "common")
diff_ex_tagwise <- exactTest(filtered_list, pair = c("No Dementia" , "Dementia"), dispersion = "tagwise")
diff_ex_poisson <- exactTest(filtered_list, pair = c("No Dementia" , "Dementia"), dispersion = 1e-06)
diff_ex_auto <- exactTest(filtered_list, pair = c("No Dementia" , "Dementia"), dispersion = "auto")

names(diff_ex_auto)
diff_ex_auto$comparison # which groups have been compared
head(diff_ex_auto$table) # results table in order of your count matrix.

# store the tables for later
results_common <- diff_ex_common$table
results_tagwise <- diff_ex_tagwise$table
results_poisson <- diff_ex_poisson$table
results_auto <- diff_ex_auto$table

# sort by PValue
results_common <- results_common[order(results_common$PValue),]
results_tagwise <- results_tagwise[order(results_tagwise$PValue),]
results_poisson <- results_poisson[order(results_poisson$PValue),]
results_auto <- results_auto[order(results_auto$PValue),]

'
Get differentially expressed genes
'
# Names of genes where the p-value is significant at the 0.05 level
diff_genes_common <- rownames(results_common)[results_common$PValue <= 0.05]
diff_genes_tagwise <- rownames(results_tagwise)[results_tagwise$PValue <= 0.05]
diff_genes_poisson <- rownames(results_poisson)[results_poisson$PValue <= 0.05]
diff_genes_auto <- rownames(results_auto)[results_auto$PValue <= 0.05]

# Total number of significantly different genes
length(diff_genes_common)        # 13,360
length(diff_genes_tagwise)       # 15,323
length(diff_genes_poisson)       # 22,615
length(diff_genes_auto)          # 15,323

# Looks like "auto" picks the tagwise dispersion as expected; drop "tagwise" from here

# Here are the percentage of total number genes that are significantly different 
length(diff_genes_common)/nrow(results_common) * 100        # 47.8%
length(diff_genes_auto)/nrow(results_auto) * 100            # 54.9%
length(diff_genes_poisson)/nrow(results_poisson) * 100      # 81.0%

# decideTestsDGE classifies the fold change metrics in the table as up, down, or NS;
# uses correction for repeated measures that controls the false discovery rate 
# (expected proportion of false discoveries among rejected hypotheses).
summary(decideTestsDGE(diff_ex_common, p.value = 0.05))
summary(decideTestsDGE(diff_ex_auto, p.value = 0.05))
summary(decideTestsDGE(diff_ex_poisson, p.value = 0.05))

'
Visualize results!
'
# histograms for the top 100 differentially expressed genes
par( mfrow=c(3 ,1) )

hist(results_poisson[diff_genes_poisson[1:100],"logCPM"], 
     breaks=25, 
     xlab="Log Concentration", 
     col="darkmagenta", 
     freq=FALSE, 
     main="Poisson: Top 100")

hist(results_common[diff_genes_common[1:100],"logCPM"], 
     breaks=25, 
     xlab="Log Concentration", 
     col="darkgreen", 
     freq=FALSE, 
     main="Common: Top 100")

hist(results_auto[diff_genes_auto[1:100],"logCPM"], 
     breaks=25, 
     xlab="Log Concentration", 
     col="dodgerblue3", 
     freq=FALSE, 
     main="Tagwise: Top 100")

par(mfrow=c(1,1))

# smear plots
par( mfrow=c(3,1) )

plotSmear(filtered_list, 
          de.tags=diff_genes_poisson, 
          main="Poisson", 
          pair = c("No Dementia" , "Dementia"),
          cex = .35,
          xlab="Log Concentration", 
          ylab="Log Fold-Change")
abline(h = c(-2, 2), col = "darkmagenta", lwd = 3)

plotSmear(filtered_list, 
          de.tags=diff_genes_common, 
          main="Common",
          pair = c("No Dementia" , "Dementia"),
          cex = .35,
          xlab="Log Concentration", 
          ylab="Log Fold-Change")
abline(h=c(-2,2), col="darkmagenta", lwd = 3)

plotSmear(filtered_list, 
          de.tags=diff_genes_auto, 
          main="Tagwise",
          pair = c("No Dementia" , "Dementia"),
          cex = .35,
          xlab="Log Concentration", 
          ylab="Log Fold-Change")
abline(h=c(-2,2), col="darkmagenta", lwd = 3)

par(mfrow=c(1,1))

# heatmap

# get top differentially expressed tags and subset by the diff_genes list
diff_tags_auto <- topTags(diff_ex_auto, n=2000000)

sig_genes <- rownames(diff_tags_auto)[(abs(diff_tags_auto$table$logFC > 2.)) & (diff_tags_auto$table$FDR < 0.05)]

heatmap.2(data.matrix(DGE_list[sig_genes,]))