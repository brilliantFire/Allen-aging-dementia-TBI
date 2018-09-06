##
'
Neuropathological Data Preparation ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
05 Sep 2018 - Created script

This script loads neuropathological data and performs transformations
or replaces missing values as necessary.

See http://aging.brain-map.org/overview/home for more details about the data.

R 3.4.1
'
##

setwd('..')

'
Install packages & load libraries
'
library(data.table)

'
Load data
'
# Luminex protein, immunohistochemistry, and isoprostane quants
neuropath_data <- data.frame(fread('http://aging.brain-map.org/api/v2/data/query.csv?criteria=model::ApiTbiDonorMetric,rma::options[num_rows$eqall]'))

# sample_info
sample_info <- readRDS(file='data/sample_info.Rds')

###
'
TEMP SECTION - Load and look at clustering results
'
### HIP
hip_results01 <- readRDS(file='data/hip_clustering_results01.Rds')
hip_conditions01 <- readRDS(file='data/hip_clustering_conditions01.Rds')

hip_results02 <- readRDS(file='data/hip_clustering_results02.Rds')
hip_conditions02 <- readRDS(file='data/hip_clustering_conditions02.Rds')

### FWM
fwm_results01 <- readRDS(file='data/fwm_clustering_results01.Rds')
fwm_conditions01 <- readRDS(file='data/fwm_clustering_conditions01.Rds')

fwm_results02 <- readRDS(file='data/fwm_clustering_results02.Rds')
fwm_conditions02 <- readRDS(file='data/fwm_clustering_conditions02.Rds')
