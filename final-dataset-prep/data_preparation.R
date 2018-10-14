##
'
Final Dataset Preparation ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
11 Sep 2018 - Script created
08 Oct 2018 - Changed N/A category in apo_4e_allele variable
11 Oct 2018 - Finalized script
12 Oct 2018 - Added pretty plots for figures
13 Oct 2018 - Correlograms

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
library(plyr)         # more of the tidyverse!
library(dplyr)        # select_if
library(matrixStats)  # colMedians
library(Amelia)       # missmap
library(mice)         # predictive mean matching
library(ggplot2)      # pretty plots for figures
library(corrplot)     # correlograms
library(RColorBrewer) # colors!

'
Load data, missmap, & drop variables
'
data <- read.csv('data/final_dataset_unprepped_CSV.csv')

# make a copy of original
data_copy <- data

# missmap for predictors
drops_missmap <- c('X', 'obs_weight', 'act_demented')
data_missmap <- data[ , !names(data) %in% drops_missmap]

png('data/missmap_allvariables_v.2.png', height=500, width=800)
missmap(data_missmap, col=c('gray80', 'royalblue4'), main='', x.cex=0.3, y.cex=1, margins=c(6.5,2.0))
dev.off()

# drop missing, zero, severely unbalanced, and duplicate variables
drops <- c('isoprostane_pg_per_mg.FWM',
           'isoprostane_pg_per_mg.FWM',
           'il_1b_pg_per_mg.FWM',
           'il_1b_pg_per_mg.TCx',
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

# make a copy of the dataset with dropped variables
data_drops <- data

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

# selected plots with and without imputed values
par(mfrow=c(2,2)) 
# ptau_ng_per_mg.PCx
plot(data$ptau_ng_per_mg.PCx, pch=22, col='royalblue4', bg='royalblue3', ylab='ptau_ng_per_mg.PCx')
points(data_mice$imp$ptau_ng_per_mg.PCx, pch=24, col='maroon4', bg='maroon3')
legend('topleft', c('observed', 'imputed'), pch=c(15,17), col=c('royalblue3', 'maroon3'))

# isoprostane_ng_per_mg.TCx
plot(data$isoprostane_pg_per_mg.TCx, pch=22, col='royalblue4', bg='royalblue3', ylab='isoprostane_ng_per_mg.TCx')
points(data_mice$imp$isoprostane_pg_per_mg.TCx, pch=24, col='maroon4', bg='maroon3')

# a_syn_pg_per_mg.HIP
plot(data$a_syn_pg_per_mg.HIP, pch=22, col='royalblue4', bg='royalblue3', ylab='a_syn_pg_per_mg.HIP')
points(data_mice$imp$a_syn_pg_per_mg.HIP, pch=24, col='maroon4', bg='maroon3')

# vegf_pg_per_mg.FWM
plot(data$vegf_pg_per_mg.FWM, pch=22, col='royalblue4', bg='royalblue3', ylab='vegf_pg_per_mg.FWM')
points(data_mice$imp$vegf_pg_per_mg.FWM, pch=24, col='maroon4', bg='maroon3')

# add back weights & target
data_imp$obs_weight <- obs_weight
data_imp$act_demented <- act_demented

write.csv(data_imp, file='data/final_dataset_cart_imputed.csv')
saveRDS(data_imp, file='data/final_dataset_cart_imputed.Rds')

'
Pretty plots for figures
'
# at8 in fft versus ffpe samples
at8_fft_hip_means <- ddply(data_copy, 'act_demented', summarise, grp.mean=mean(ihc_at8.HIP, na.rm=TRUE))
at8_ffpe_hip_means <- ddply(data_copy, 'act_demented', summarise, grp.mean=mean(ihc_at8_ffpe.HIP, na.rm=TRUE))

at8_fft_fwm_means <- ddply(data_copy, 'act_demented', summarise, grp.mean=mean(ihc_at8.FWM, na.rm=TRUE))
at8_ffpe_fwm_means <- ddply(data_copy, 'act_demented', summarise, grp.mean=mean(ihc_at8_ffpe.FWM, na.rm=TRUE))

at8_fft_pcx_means <- ddply(data_copy, 'act_demented', summarise, grp.mean=mean(ihc_at8.PCx, na.rm=TRUE))
at8_ffpe_pcx_means <- ddply(data_copy, 'act_demented', summarise, grp.mean=mean(ihc_at8_ffpe.PCx, na.rm=TRUE))

at8_fft_tcx_means <- ddply(data_copy, 'act_demented', summarise, grp.mean=mean(ihc_at8.TCx, na.rm=TRUE))
at8_ffpe_tcx_means <- ddply(data_copy, 'act_demented', summarise, grp.mean=mean(ihc_at8_ffpe.TCx, na.rm=TRUE))

# TCx example fft
png('data/tcx_at8_fft_densities_orig.png', height=500, width=700)
ggplot(data_copy, aes(x=ihc_at8.TCx, fill=act_demented)) +
  geom_density(alpha=0.3) +
  geom_vline(data=at8_fft_tcx_means, aes(xintercept=grp.mean, color=act_demented),
                  linetype='dashed', size=1, show.legend = F) +
  labs(x='% of Area', y='density') +
  guides(fill=guide_legend(title='')) +
  theme(axis.title=element_text(size=24), axis.text=element_text(size=22), legend.text=element_text(size=22))
dev.off()

# TCx exampleffpe
png('data/tcx_at8_ffpe_densities_orig.png', height=500, width=700)
ggplot(data_copy, aes(x=ihc_at8_ffpe.TCx, fill=act_demented)) +
  geom_density(alpha=0.3) +
  geom_vline(data=at8_ffpe_tcx_means, aes(xintercept=grp.mean, color=act_demented),
                  linetype='dashed', size=1, show.legend = F) +
  labs(x='% of Area', y='density') +
  guides(fill=guide_legend(title='')) +
  theme(axis.title=element_text(size=24), axis.text=element_text(size=22), legend.text=element_text(size=22))
dev.off()

# HIP ptau/tau ratio and components
ptau_hip_means <- ddply(data_copy, 'act_demented', summarise, grp.mean=mean(ptau_ng_per_mg.HIP, na.rm=TRUE))
tau_hip_means <- ddply(data_copy, 'act_demented', summarise, grp.mean=mean(tau_ng_per_mg.HIP, na.rm=TRUE))
tau_ratio_hip_means <- ddply(data_copy, 'act_demented', summarise, grp.mean=mean(ptau_over_tau_ratio.HIP, na.rm=TRUE))

png('data/hip_ptau_densities_orig.png', height=500, width=700)
ggplot(data_copy, aes(x=ptau_ng_per_mg.HIP, fill=act_demented)) +
  geom_density(alpha=0.3) +
  geom_vline(data=ptau_hip_means, aes(xintercept=grp.mean, color=act_demented),
                  linetype='dashed', size=1, show.legend = F) +
  labs(x='ptau (ng/mg)', y='density') +
  guides(fill=guide_legend(title='')) +
  theme(axis.title=element_text(size=24), axis.text=element_text(size=22), legend.text=element_text(size=22))
dev.off()

png('data/hip_tau_densities_orig.png', height=500, width=700)
ggplot(data_copy, aes(x=tau_ng_per_mg.HIP, fill=act_demented)) +
  geom_density(alpha=0.3) +
  geom_vline(data=tau_hip_means, aes(xintercept=grp.mean, color=act_demented),
                  linetype='dashed', size=1, show.legend = F) +
  labs(x='tau (ng/mg)', y='density') +
  guides(fill=guide_legend(title='')) +
  theme(axis.title=element_text(size=24), axis.text=element_text(size=22), legend.text=element_text(size=22))
dev.off()

png('data/hip_tau_ratio_densities_orig.png', height=500, width=700)
ggplot(data_copy, aes(x=ptau_over_tau_ratio.HIP, fill=act_demented)) +
  geom_density(alpha=0.3) +
  geom_vline(data=tau_ratio_hip_means, aes(xintercept=grp.mean, color=act_demented),
                  linetype='dashed', size=1, show.legend = F) +
  labs(x='ptau/tau ratio', y='density') +
  guides(fill=guide_legend(title='')) +
  theme(axis.title=element_text(size=24), axis.text=element_text(size=22), legend.text=element_text(size=22))
dev.off()

'
Correlograms
'
# subset by brain region
hip <- c('ihc_a_syn.HIP',
         'ihc_at8_ffpe.HIP',
         'ihc_tau2_ffpe.HIP',
         'ihc_a_beta_ffpe.HIP',
         'ihc_ptdp_43_ffpe.HIP',
         'ihc_iba1_ffpe.HIP',
         'ihc_gfap_ffpe.HIP',
         'ptau_ng_per_mg.HIP',
         'vegf_pg_per_mg.HIP',
         'tnf_a_pg_per_mg.HIP',
         'tau_ng_per_mg.HIP',
         'il_10_pg_per_mg.HIP',
         'il_6_pg_per_mg.HIP',
         'il_1b_pg_per_mg.HIP',
         'rantes_pg_per_mg.HIP',
         'ab40_pg_per_mg.HIP',
         'a_syn_pg_per_mg.HIP',
         'ifn_g_pg_per_mg.HIP',
         'mcp_1_pg_per_mg.HIP',
         'bdnf_pg_per_mg.HIP',
         'mip_1a_pg_per_mg.HIP',
         'il_7_pg_per_mg.HIP',
         'ab42_pg_per_mg.HIP',
         'gene_cluster01.HIP',
         'gene_cluster02.HIP',
         'gene_cluster03.HIP',
         'education_years',
         'age_at_first_tbi',
         'cerad',
         'num_tbi_w_loc',
         'braak',
         'nia_reagan')

fwm <- c('ihc_a_syn.FWM',
         'ihc_at8_ffpe.FWM',
         'ihc_tau2_ffpe.FWM',
         'ihc_a_beta_ffpe.FWM',
         'ihc_ptdp_43_ffpe.FWM',
         'ihc_iba1_ffpe.FWM',
         'ihc_gfap_ffpe.FWM',
         'ptau_ng_per_mg.FWM',
         'vegf_pg_per_mg.FWM',
         'tnf_a_pg_per_mg.FWM',
         'tau_ng_per_mg.FWM',
         'il_10_pg_per_mg.FWM',
         'il_6_pg_per_mg.FWM',
         'il_1b_pg_per_mg.FWM',
         'rantes_pg_per_mg.FWM',
         'ab40_pg_per_mg.FWM',
         'a_syn_pg_per_mg.FWM',
         'ifn_g_pg_per_mg.FWM',
         'mcp_1_pg_per_mg.FWM',
         'bdnf_pg_per_mg.FWM',
         'mip_1a_pg_per_mg.FWM',
         'il_7_pg_per_mg.FWM',
         'ab42_pg_per_mg.FWM',
         'gene_cluster01.FWM',
         'gene_cluster02.FWM',
         'education_years',
         'age_at_first_tbi',
         'cerad',
         'num_tbi_w_loc',
         'braak',
         'nia_reagan')

pcx <- c('ihc_a_syn.PCx',
         'ihc_at8_ffpe.PCx',
         'ihc_tau2_ffpe.PCx',
         'ihc_a_beta_ffpe.PCx',
         'ihc_ptdp_43_ffpe.PCx',
         'ihc_iba1_ffpe.PCx',
         'ihc_gfap_ffpe.PCx',
         'ptau_ng_per_mg.PCx',
         'vegf_pg_per_mg.PCx',
         'tnf_a_pg_per_mg.PCx',
         'tau_ng_per_mg.PCx',
         'il_10_pg_per_mg.PCx',
         'isoprostane_pg_per_mg.PCx',
         'il_6_pg_per_mg.PCx',
         'il_1b_pg_per_mg.PCx',
         'rantes_pg_per_mg.PCx',
         'ab40_pg_per_mg.PCx',
         'a_syn_pg_per_mg.PCx',
         'ifn_g_pg_per_mg.PCx',
         'mcp_1_pg_per_mg.PCx',
         'bdnf_pg_per_mg.PCx',
         'mip_1a_pg_per_mg.PCx',
         'il_7_pg_per_mg.PCx',
         'ab42_pg_per_mg.PCx',
         'gene_cluster01.PCx',
         'gene_cluster02.PCx',
         'gene_cluster03.PCx',
         'education_years',
         'age_at_first_tbi',
         'cerad',
         'num_tbi_w_loc',
         'braak',
         'nia_reagan')

tcx <- c('ihc_a_syn.TCx',
         'ihc_at8_ffpe.TCx',
         'ihc_tau2_ffpe.TCx',
         'ihc_a_beta_ffpe.TCx',
         'ihc_ptdp_43_ffpe.TCx',
         'ihc_iba1_ffpe.TCx',
         'ihc_gfap_ffpe.TCx',
         'ptau_ng_per_mg.TCx',
         'vegf_pg_per_mg.TCx',
         'tnf_a_pg_per_mg.TCx',
         'tau_ng_per_mg.TCx',
         'il_10_pg_per_mg.TCx',
         'isoprostane_pg_per_mg.TCx',
         'il_6_pg_per_mg.TCx',
         'il_1b_pg_per_mg.TCx',
         'rantes_pg_per_mg.TCx',
         'ab40_pg_per_mg.TCx',
         'a_syn_pg_per_mg.TCx',
         'ifn_g_pg_per_mg.TCx',
         'mcp_1_pg_per_mg.TCx',
         'bdnf_pg_per_mg.TCx',
         'mip_1a_pg_per_mg.TCx',
         'il_7_pg_per_mg.TCx',
         'ab42_pg_per_mg.TCx',
         'gene_cluster01.TCx',
         'gene_cluster02.TCx',
         'education_years',
         'age_at_first_tbi',
         'cerad',
         'num_tbi_w_loc',
         'braak',
         'nia_reagan')

all_vars <- data_drops[, unlist(lapply(data_drops, is.numeric))]

hip_vars <- data_drops[ , names(data_drops) %in% hip]
fwm_vars <- data_drops[ , names(data_drops) %in% fwm]
pcx_vars <- data_drops[ , names(data_drops) %in% pcx]
tcx_vars <- data_drops[ , names(data_drops) %in% tcx]

# HIP
png('data/hip_pearson_correlogram.png', height=900, width=900)
corrplot(cor(hip_vars, use = 'complete.obs'), tl.cex=1.3, tl.col='black', order='hclust',
         col=brewer.pal(n=10, name='PiYG'))
dev.off()

# FWM
png('data/fwm_pearson_correlogram.png', height=900, width=900)
corrplot(cor(fwm_vars, use = 'complete.obs'), tl.cex=1.3, tl.col='black', order='hclust',
         col=brewer.pal(n=10, name='PiYG'))
dev.off()

# PCx
png('data/pcx_pearson_correlogram.png', height=900, width=900)
corrplot(cor(pcx_vars, use = 'complete.obs'), tl.cex=1.3, tl.col='black', order='hclust',
         col=brewer.pal(n=10, name='PiYG'))
dev.off()

# TCx
png('data/tcx_pearson_correlogram.png', height=900, width=900)
corrplot(cor(tcx_vars, use = 'complete.obs'), tl.cex=1.3, tl.col='black', order='hclust',
         col=brewer.pal(n=10, name='PiYG'))
dev.off()

# All
png('data/all_pearson_correlogram.png', height=1000, width=1000)
corrplot(cor(all_vars, use='complete.obs'), tl.cex=0.7, tl.col='black', order='hclust',
         col=brewer.pal(n=10, name='PiYG'))
dev.off()