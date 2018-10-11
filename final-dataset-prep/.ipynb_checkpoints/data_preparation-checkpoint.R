##
'
Final Dataset Preparation ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
11 Sep 2018 - Script created
08 Oct 2018 - Changed N/A category in apo_4e_allele variable
11 Oct 2018 - Finalized script

This script prepares the dataset constructed in `feature_engineering_dataset_construction.R`
for modeling, including missing value imputation by classification & regression trees (CART)
using the mice package.

See http://aging.brain-map.org/overview/home for more details about the data.

R 3.5.1
'
##

setwd('..')
set.seed(555)

'
Load libraries
'
library(dplyr)        # select_if
library(matrixStats)  # colMedians
library(mice)         # predictive mean matching

'
Load data, drop variables
'
data <- read.csv('data/final_dataset_unprepped_CSV.csv')

# drop missing, zero, severely unbalanced, and duplicate variables
drops <- c('isoprostane_pg_per_mg.HIP',
           'isoprostane_pg_per_mg.FWM',
           'il_1b_pg_per_mg.FWM',
           'il_1b_pg_per_mg.PCx',
           'il_1b_pg_per_mg.TCx',
           'il_4_pg_per_mg.HIP',
           'ab42_over_ab40_ratio.HIP',
           'ab42_over_ab40_ratio.FWM',
           'ab42_over_ab40_ratio.PCx',
           'ab42_over_ab40_ratio.TCx',
           'ptau_over_tau_ratio.HIP',
           'ptau_over_tau_ratio.FWM',
           'ptau_over_tau_ratio.PCx',
           'ptau_over_tau_ratio.TCx',
           'ihc_at8.HIP',
           'ihc_at8.FWM',
           'ihc_at8.PCx',
           'ihc_at8.TCx',
           'ihc_a_beta.HIP',
           'ihc_a_beta.FWM',
           'ihc_a_beta.PCx',
           'ihc_a_beta.TCx',
           'race',
           'hispanic',
           'X')

data <- data[ , !names(data) %in% drops]

# structure list
str(data, list.len=ncol(data))

'
Mean/Median Comparisons for Numeric Variables
'
# select all numeric variables
numerics <- select_if(data, is.double)

# means, medians, & % difference of the median vs the mean
variables <- names(numerics)
means <- colMeans(numerics, na.rm=TRUE)
medians <- colMedians(as.matrix(numerics), na.rm=TRUE)
per_diff <- round(((means - medians)/means)*100, 2)

# combine
numeric_metrics <- data.frame(cbind(means, medians, per_diff))
colnames(numeric_metrics) <- c('Means', 'Medians', 'Percent Diff')
rownames(numeric_metrics) <- variables

'
Replace age Ranges with centers
'
data$age <- as.character(data$age)

data$age[data$age == '100+'] <- 100
data$age[data$age == '90-94'] <- 92
data$age[data$age == '95-99'] <- 97

data$age <- as.numeric(data$age)

'
Reorder Levels for Some Categorical Variables
'
# longest_loc_duration
loc_durations <- c('Unknown or N/A',
                   '< 10 sec',
                   '10 sec - 1 min',
                   '1-2 min',
                   '3-5 min',
                   '6-9 min',
                   '10 min - 1 hr',
                   '> 1 hr')

data$longest_loc_duration <- factor(data$longest_loc_duration, 
                                    ordered=TRUE,
                                    levels=loc_durations)

# dsm_iv_clinical_diagnosis
dsm_diagnoses <- c('No Dementia',
                   "Alzheimer's Disease Type",
                   'Vascular',
                   'Other or Unknown Cause',
                   'Other Medical',
                   'Multiple Etiologies')

data$dsm_iv_clinical_diagnosis <- factor(data$dsm_iv_clinical_diagnosis, levels=dsm_diagnoses)

# nincds_arda_diagnosis
nincds_diagnoses <- c('No Dementia',
                      "Possible Alzheimer's Disease",
                      "Probable Alzheimer's Disease",
                      'Dementia, Type Unknown')

data$nincds_arda_diagnosis <- factor(data$nincds_arda_diagnosis, levels=nincds_diagnoses)

'
Factorize Other Categorical Variables (including target, act_demented)
'
# sex
data$sex <- factor(data$sex)

# apo_e4_allele (first replace N/A with most frequent category 'N')
data$apo_e4_allele[data$apo_e4_allele == 'N/A'] <- 'N'
data$apo_e4_allele <- factor(data$apo_e4_allele)

# ever_tbi_w_loc
data$ever_tbi_w_loc <- factor(data$ever_tbi_w_loc)

# cerad
data$cerad <- factor(data$cerad)

# braak
data$braak <- factor(data$braak)

# nia_reagan
data$nia_reagan <- factor(data$nia_reagan)

# act_demented
data$act_demented <- factor(data$act_demented)

'
Imputation of Missing Values by Classification & Regression Trees (CART)
'
# Remove observation weights
obs_weight <- data$obs_weight
data$obs_weight <- NULL

# Remove act_demented
act_demented <- data$act_demented
data$act_demented <- NULL

# Impute missing values
data_mice <- mice(data, m=1, method='cart')
data_imp <- complete(data_mice)

# add back weights & target
data_imp$obs_weight <- obs_weight
data_imp$act_demented <- act_demented

write.csv(data_imp, file='data/final_dataset_cart_imputed.csv')
saveRDS(data_imp, file='data/final_dataset_cart_imputed.Rds')
