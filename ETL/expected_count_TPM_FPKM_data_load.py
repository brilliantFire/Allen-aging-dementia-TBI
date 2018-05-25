##
'''
Load Unnormalized Expected Count, TPM, & FPKM Data ~ Python version
Rebecca Vislay Wade
01 May 2018
25 May 2018 - updated i/o to work with GitHub repo

Loads gene expression data for 377 brain samples from the *Aging, Dementia, & TBI
Study* from the Allen Institute for Brain Science. Saves resulting dataframes to
.csv files.

See http://aging.brain-map.org/overview/home for more details about the data.

Python 3.6.5
'''
##

# Some imports
import pandas as pd
import numpy as np
import time
import os

# load TBI_data_files.csv
data_files = pd.read_csv('http://aging.brain-map.org/data/tbi_data_files.csv')

# Series of links to files containing gene level estimated count/TPM/FPKM for each sample 
data_links = data_files['gene_level_fpkm_file_link']

# Grab RNA-seq sample IDs (these will be the column names in the final DFs)
sample_IDs = data_files['rnaseq_profile_id']

# Initialize DFs
TPM = pd.DataFrame()
FPKM = pd.DataFrame()
raw_read_counts = pd.DataFrame()

# start a timer for dataframe construction
start = time.clock()

# loop through list of links
for sample in range(0, len(data_links)):
    print('Now adding Sample #', sample)
    # get data file for the sample
    url = 'http://aging.brain-map.org' + str(data_links[sample])
    sample_data = pd.read_table(url)
    # get gene_id from the first sample file (they are all the same)
    if sample == 0:
        TPM['gene_id'] = sample_data['gene_id']
        FPKM['gene_id'] = sample_data['gene_id']
        raw_read_counts['gene_id'] = sample_data['gene_id']
    # add TPM, FPKM, and expected_count to respective DataFrames
    TPM[str(sample_IDs[sample])] = sample_data['TPM']
    FPKM[str(sample_IDs[sample])] = sample_data['FPKM']
    raw_read_counts[str(sample_IDs[sample])] = sample_data['expected_count']

# set gene_id as indices
TPM = TPM.set_index('gene_id')
FPKM = FPKM.set_index('gene_id')
raw_read_counts = raw_read_counts.set_index('gene_id')

# stop timer
stop = time.clock()
duration = stop - start

# save to csv for later
os.chdir('..')
TPM.to_csv('data/TPM.csv')
FPKM.to_csv('data/FPKM.csv')
raw_read_counts.to_csv('data/raw_read_counts.csv')

print('All data loaded & saved! Loop duration:', np.round(duration, 2), 'seconds')