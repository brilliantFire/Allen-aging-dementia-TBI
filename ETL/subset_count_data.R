##
'
Differential Gene Expression Analysis ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
24 May 2018
25 May 2018 - updated i/o to work with GitHub repo

This script splits the raw RNA-seq count data (loaded using the 
"expected_count_TPM_FPKM_data_load" script) into groups both by brain
region and by donor sex. Outputs CSVs for each group that have the count
data plus dementia status in the last row ([50284,]).

See http://aging.brain-map.org/overview/home for more details about the data.

R 3.4.1
'
##

'
Load libraries
'
library(data.table)    # i/o

'
Load data
'
setwd('..')
raw_read_counts <- data.frame(read.csv('data/raw_read_counts.csv'))

# make 'gene_id' row names
rownames(raw_read_counts) <- raw_read_counts$gene_id
raw_read_counts$gene_id <- NULL

# Remove X from column names
colnames(raw_read_counts) <- substring(colnames(raw_read_counts), 2)

# Get donor info
donor_files <- data.frame(fread('http://aging.brain-map.org/api/v2/data/query.csv?criteria=model::ApiTbiDonorDetail,rma::options[num_rows$eqall]'))

# Get TBI_data_files.csv
data_files <- data.frame(fread('http://aging.brain-map.org/data/tbi_data_files.csv'))

# merge on donor_id
sample_info <- merge(data_files, donor_files, by = 'donor_id', all = TRUE)

'
Makes lists of RNA-seq profile IDs for each brain region & sex
'
# get list of samples in each brain region group
hip_list <- sample_info$rnaseq_profile_id[which(sample_info$structure_acronym == 'HIP')]
fwm_list <- sample_info$rnaseq_profile_id[which(sample_info$structure_acronym == 'FWM')]
pcx_list <- sample_info$rnaseq_profile_id[which(sample_info$structure_acronym == 'PCx')]
tcx_list <- sample_info$rnaseq_profile_id[which(sample_info$structure_acronym == 'TCx')]

male_list <- sample_info$rnaseq_profile_id[which(sample_info$sex == 'M')]
female_list <- sample_info$rnaseq_profile_id[which(sample_info$sex == 'F')]

'
Loop through lists and collect samples in individual matrices
'
# initialize matrices
hip_counts <- matrix()
fwm_counts <- matrix()
pcx_counts <- matrix()
tcx_counts <- matrix()

male_counts <- matrix()
female_counts <- matrix()

# loops to collect samples
for (i in 1:length(hip_list)){
    sample <- raw_read_counts[which(colnames(raw_read_counts) == hip_list[i])] 
    hip_counts <- cbind(hip_counts, sample)
}

for (i in 1:length(fwm_list)){
    sample <- raw_read_counts[which(colnames(raw_read_counts) == fwm_list[i])] 
    fwm_counts <- cbind(fwm_counts, sample)
}

for (i in 1:length(pcx_list)){
    sample <- raw_read_counts[which(colnames(raw_read_counts) == pcx_list[i])] 
    pcx_counts <- cbind(pcx_counts, sample)
}

for (i in 1:length(tcx_list)){
    sample <- raw_read_counts[which(colnames(raw_read_counts) == tcx_list[i])] 
    tcx_counts <- cbind(tcx_counts, sample)
}

for (i in 1:length(male_list)){
    sample <- raw_read_counts[which(colnames(raw_read_counts) == male_list[i])] 
    male_counts <- cbind(male_counts, sample)
}

for (i in 1:length(female_list)){
    sample <- raw_read_counts[which(colnames(raw_read_counts) == female_list[i])] 
    female_counts <- cbind(female_counts, sample)
}

# remove extra columns
hip_counts <- hip_counts[,-1]
fwm_counts <- fwm_counts[,-1]
pcx_counts <- pcx_counts[,-1]
tcx_counts <- tcx_counts[,-1]

male_counts <- male_counts[,-1]
female_counts <- female_counts[,-1]

'
Get dementia statuses 
'
hip_status <- sample_info$act_demented[which(sample_info$structure_acronym == 'HIP')]
fwm_status <- sample_info$act_demented[which(sample_info$structure_acronym == 'FWM')]
pcx_status <- sample_info$act_demented[which(sample_info$structure_acronym == 'PCx')]
tcx_status <- sample_info$act_demented[which(sample_info$structure_acronym == 'TCx')]

male_status <- sample_info$act_demented[which(sample_info$sex == 'M')]
female_status <- sample_info$act_demented[which(sample_info$sex == 'F')]

'
Construct & save data frames
'
# Construct dataframes with counts & dementia statuses
hip_counts <- rbind(hip_counts, hip_status)
fwm_counts <- rbind(fwm_counts, fwm_status)
pcx_counts <- rbind(pcx_counts, pcx_status)
tcx_counts <- rbind(tcx_counts, tcx_status)

male_counts <- rbind(male_counts, male_status)
female_counts <- rbind(female_counts, female_status)

# write to csv in ../data
write.csv(hip_counts, 'data/hip_counts_status.csv')
write.csv(fwm_counts, 'data/fwm_counts_status.csv')
write.csv(pcx_counts, 'data/pcx_counts_status.csv')
write.csv(tcx_counts, 'data/tcx_counts_status.csv')

write.csv(male_counts, 'data/male_counts_status.csv')
write.csv(female_counts, 'data/female_counts_status.csv')
