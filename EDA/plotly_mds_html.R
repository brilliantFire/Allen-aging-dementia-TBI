##
'
Differential Gene Expression Analysis ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
23 May 2018

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

# load counts from local copy
raw_read_counts <- data.frame(read.csv('C:/Users/Rebecca/Documents/NU_MSPA/thesis/data/raw_read_counts.csv'))

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

# plotly shaded by sex
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
       title = 'MDS Plot of Gene Expression Profiles, Shaded by Donor Sex')
htmlwidgets::saveWidget(as_widget(p_sex), "by_sex.html")

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
       title = 'MDS Plotof Gene Expression Profiles, Shaded by Brain Region')
htmlwidgets::saveWidget(as_widget(p_region), "by_region.html")
