##
'
Load Unnormalized Expected Count, TPM, & FPKM Data ~ R version
Rebecca Vislay Wade
01 May 2018

Loads gene expression data for 377 brain samples from the *Aging, Dementia, & TBI
Study* from the Allen Institute for Brain Science. Saves resulting dataframes to
.csv files.

See http://aging.brain-map.org/overview/home for more details about the data.

R 3.4.2
'
##

# load libraries
library(data.table)

# Get TBI_data_files.csv
data_files <- data.frame(fread('http://aging.brain-map.org/data/tbi_data_files.csv'))

# Series of links to normalized FMPK files
data_links <- sample_info['gene_level_fpkm_file_link']

# Grab the list of sample IDs
sample_IDs <- sample_info['rnaseq_profile_id']

# Initialize dataframes
TPM <- data.frame(matrix(nrow = 50283))
FPKM <- data.frame(matrix(nrow = 50283))
raw_read_counts <- data.frame(matrix(nrow = 50283))

# start a timer for DF construction
start <- Sys.time()

# loops through sample files and adds raw read counts ('expected_counts') for each
# sample to a dataframe with the 'rnaseq_profile_id' as header
for (sample in 1:nrow(data_links)){
    print(paste('Sample', sample, 'being added.'))
    flush.console()     # Forces print in loop
    url <- paste('http://aging.brain-map.org', toString(data_links[sample,1]),sep='')
    sample_data <- data.frame(fread(url))
    # Add 'gene_id' from the first file to each dataframe (it's the same for each file)
    if (sample == 1){
        TPM['gene_id'] <- sample_data['gene_id']
        FPKM['gene_id'] <- sample_data['gene_id']
        raw_read_counts['gene_id'] <- sample_data['gene_id']
        }
    TPM[toString(sample_IDs[sample,1])] <- sample_data['TPM']
    FPKM[toString(sample_IDs[sample,1])] <- sample_data['FPKM']
    raw_read_counts[toString(sample_IDs[sample,1])] <- sample_data['expected_count']
    }

# stop timer, calculate running time
stop <- Sys.time()
duration <- stop-start

# drop extra column
TPM <- TPM[,-1]
FPKM <- FPKM[,-1]
raw_read_counts <- raw_read_counts[,-1]

# Write to csv for later
write.csv(raw_read_counts,'raw_read_counts.csv')
write.csv(TPM,'TPM.csv')
write.csv(FPKM,'FPKM.csv')

print(paste('All samples loaded & saved! Loop duration:', round(duration, 2), 'minutes'))