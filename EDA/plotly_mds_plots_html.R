##
'
Differential Gene Expression Analysis ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
23 May 2018
25 May 2018 - updated i/o to work with GitHub repo

This script loads, cleans, filters, and normalizes counts. 3D plotly
multidimensional scaling plots are created and saved as html.

See http://aging.brain-map.org/overview/home for more details about the data.

R 3.4.1
'
##
# load libraries
library(data.table)    # i/o
library(edgeR)         # DGE (also loads dependencies)
library(plotly)        # Awesome 3D interactive plots 

# load counts from copy in ../data
setwd('..')
raw_read_counts <- data.frame(read.csv('data/raw_read_counts.csv'))

# make 'gene_id' row names
rownames(raw_read_counts) <- raw_read_counts$gene_id
raw_read_counts$gene_id <- NULL

# Remove X from column names
colnames(raw_read_counts) <- substring(colnames(raw_read_counts), 2)

# get donor info
donor_files <- data.frame(fread('http://aging.brain-map.org/api/v2/data/query.csv?criteria=model::ApiTbiDonorDetail,rma::options[num_rows$eqall]'))

# Get TBI_data_files.csv
data_files <- data.frame(fread('http://aging.brain-map.org/data/tbi_data_files.csv'))

# merge on donor_id
sample_info <- merge(data_files, donor_files, by = 'donor_id', all = TRUE)

# get target variable act_demented
dementia_status <- t(sample_info$act_demented)

# puts counts into 'counts' & group into 'samples'
DGE_list <- DGEList(raw_read_counts, group = dementia_status)

# filter counts
filtered_counts <- DGE_list[rowSums(1e+06*DGE_list$counts/expandAsMatrix(DGE_list$samples$lib.size, dim(raw_read_counts))>2)>=10,]

# calculate normalization factors
norm_counts <- calcNormFactors(filtered_counts, method = 'TMM')

# MDS
logFC_500 <- plotMDS(norm_counts,
                     method = "logFC",
                     top = 500,
                     ndim = 4)

# make a dataframe with 4D coordinates + brain region, sex, & dementia status from sample_info
for_plots <- data.frame(logFC_500$cmdscale.out) 
colnames(for_plots) <- c('dim1','dim2','dim3','dim4')
for_plots$brain_region <- sample_info$structure_acronym
for_plots$sex <- sample_info$sex
for_plots$dementia_status <- sample_info$act_demented

# plotly plot shaded by sex & save as html widgets in ../data
p_sex <- plot_ly(for_plots, 
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
       title = 'Gene Expression Profile MDS Plot, Shaded by Donor Sex')

# saveWidgets has problems saving to dirs other than the current one
# here's a workaround I found here:
# https://stackoverflow.com/questions/41399795/savewidget-from-htmlwidget-in-r-cannot-save-html-file-in-another-folder
temp_path <- 'data/plotly_mds_by_sex.html'
htmlwidgets::saveWidget(as_widget(p_sex), file.path(normalizePath(dirname(temp_path)),basename(temp_path)))

# by brain region
p_region <- plot_ly(for_plots, 
                    x = ~dim1, 
                    y = ~dim2, 
                    z = ~dim3, 
                    color = ~brain_region,
                    width = 700,
                    height = 700) %>%
add_markers() %>%
layout(scene = list(xaxis = list(title = 'dim-1'),
                    yaxis = list(title = 'dim-2'),
                    zaxis = list(title = 'dim-3')),
       title = 'Gene Expression Profile MDS Plot, Shaded by Brain Region')
temp_path <- 'data/plotly_MDS_by_region.html'
htmlwidgets::saveWidget(as_widget(p_sex), file.path(normalizePath(dirname(temp_path)),basename(temp_path)))