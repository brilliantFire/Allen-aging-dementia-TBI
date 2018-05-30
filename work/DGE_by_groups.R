##
'
Differential Gene Expression Analysis ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
25 May 2018

This script identifies sets of differentially expressed genes in subsets
of the Allen Aging, Dementia, and TBI RNA-seq dataset.

See http://aging.brain-map.org/overview/home for more details about the data.

R 3.4.1
'
##

## NOT WORKING YET ##

'
Load libraries
'
library(data.table)    # i/o
library(edgeR)         # DGE (also loads dependencies)
library(plotly)        # Awesome 3D interactive plots 

'
Load data
'
setwd('..')
hip_counts <- read.csv('data/hip_counts_status.csv')

# make DGEList objects
dgeList_hip <- DGEList(hip_counts, group = hip_status)
dgeList_fwm <- DGEList(fwm_counts, group = fwm_status)
dgeList_pcx <- DGEList(pcx_counts, group = pcx_status)
dgeList_tcx <- DGEList(tcx_counts, group = tcx_status)

dgeList_male <- DGEList(male_counts, group = male_status)
dgeList_female <- DGEList(female_counts, group = female_status)

'
Filter counts, keeping only genes with more than 2 counts in at least 10 samples
'
filtered_hip <- dgeList_hip[rowSums(1e+06*dgeList_hip$counts/expandAsMatrix(dgeList_hip$samples$lib.size, dim(hip_counts))>2)>=10,]
filtered_fwm <- dgeList_fwm[rowSums(1e+06*dgeList_fwm$counts/expandAsMatrix(dgeList_fwm$samples$lib.size, dim(fwm_counts))>2)>=10,]
filtered_pcx <- dgeList_pcx[rowSums(1e+06*dgeList_pcx$counts/expandAsMatrix(dgeList_pcx$samples$lib.size, dim(pcx_counts))>2)>=10,]
filtered_tcx <- dgeList_tcx[rowSums(1e+06*dgeList_tcx$counts/expandAsMatrix(dgeList_tcx$samples$lib.size, dim(tcx_counts))>2)>=10,]

filtered_male <- dgeList_male[rowSums(1e+06*dgeList_male$counts/expandAsMatrix(dgeList_male$samples$lib.size, dim(male_counts))>2)>=10,]
filtered_female <- dgeList_female[rowSums(1e+06*dgeList_female$counts/expandAsMatrix(dgeList_female$samples$lib.size, dim(female_counts))>2)>=10,]

dim(filtered_hip)       # 20311 genes
dim(filtered_fwm)       # 20646 genes
dim(filtered_pcx)       # 20597 genes
dim(filtered_tcx)       # 20587 genes

dim(filtered_male)      # 21815 genes
dim(filtered_female)    # 21357 genes

'
Between-sample normalization
'
norm_hip <- calcNormFactors(filtered_hip, method = 'TMM')
norm_fwm <- calcNormFactors(filtered_fwm, method = 'TMM')
norm_pcx <- calcNormFactors(filtered_pcx, method = 'TMM')
norm_tcx <- calcNormFactors(filtered_tcx, method = 'TMM')

norm_male <- calcNormFactors(filtered_male, method = 'TMM')
norm_female <- calcNormFactors(filtered_female, method = 'TMM')

'
MDS Plots
'
logFC_hip <- plotMDS(norm_hip, method = "logFC", top = 500, ndim = 4, main = 'hippocampus')
logFC_fwm <- plotMDS(norm_fwm, method = "logFC", top = 500, ndim = 4, main = 'forebrain white matter')
logFC_pcx <- plotMDS(norm_pcx, method = "logFC", top = 500, ndim = 4, main = 'parietal cortex')
logFC_tcx <- plotMDS(norm_tcx, method = "logFC", top = 500, ndim = 4, main = 'temporal cortex')

logFC_male <- plotMDS(norm_male, method = "logFC", top = 500, ndim = 4, main = 'males')
logFC_female <- plotMDS(norm_female, method = "logFC", top = 500, ndim = 4, main = 'females')

# make DFs for plotting
plots_hip <- data.frame(logFC_hip$cmdscale.out) 
colnames(plots_hip) <- c('dim1','dim2','dim3','dim4')
plots_hip$sex <- sample_info$sex[which(sample_info$structure_acronym == 'HIP')]
plots_hip$dementia_status <- hip_status

plots_fwm <- data.frame(logFC_fwm$cmdscale.out) 
colnames(plots_fwm) <- c('dim1','dim2','dim3','dim4')
plots_fwm$sex <- sample_info$sex[which(sample_info$structure_acronym == 'FWM')]
plots_fwm$dementia_status <- fwm_status

plots_pcx <- data.frame(logFC_pcx$cmdscale.out) 
colnames(plots_pcx) <- c('dim1','dim2','dim3','dim4')
plots_pcx$sex <- sample_info$sex[which(sample_info$structure_acronym == 'PCx')]
plots_pcx$dementia_status <- pcx_status

plots_tcx <- data.frame(logFC_tcx$cmdscale.out) 
colnames(plots_tcx) <- c('dim1','dim2','dim3','dim4')
plots_tcx$sex <- sample_info$sex[which(sample_info$structure_acronym == 'TCx')]
plots_tcx$dementia_status <- tcx_status

plots_male <- data.frame(logFC_male$cmdscale.out) 
colnames(plots_male) <- c('dim1','dim2','dim3','dim4')
plots_male$brain_region <- sample_info$structure_acronym[which(sample_info$sex == 'M')]
plots_male$dementia_status <- male_status

plots_female <- data.frame(logFC_female$cmdscale.out) 
colnames(plots_female) <- c('dim1','dim2','dim3','dim4')
plots_female$brain_region <- sample_info$structure_acronym[which(sample_info$sex == 'F')]
plots_female$dementia_status <- female_status

'
Plotly plot (look into rendering in html)
'
p_hip <- plot_ly(plots_hip, 
                 x = ~dim1, 
                 y = ~dim2, 
                 z = ~dim3, 
                 color = ~sex,
                 width = 700,
                 height = 700) %>%
add_markers() %>%
layout(scene = list(xaxis = list(title = 'dim-1'),
                    yaxis = list(title = 'dim-2'),
                    zaxis = list(title = 'dim-3')),
       title = 'Hippocampus, Shaded by Donor Sex')
htmlwidgets::saveWidget(as_widget(p_hip), "hippocampus_by_sex.html")
