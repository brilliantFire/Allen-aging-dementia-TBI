# -*- coding: utf-8 -*-
"""
~*~*~*~* Agglomerative Hierarchical Clustering Heatmap *~*~*~*~
Rebecca Vislay Wade
24 Apr 2018
--------------------------------------------------------------------------------
Uses seaborn to create a heatmap of normalized FPKM values from the Aging, 
Dementia, and TBI study by the Allen Brain Institute.
"""

# Some imports
import pandas as pd
import os
import random

from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import calinski_harabaz_score

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

# Set the working directory
#os.chdir('C:\\Users\\rlvis\\work_MSPA\\thesis') # desktop
os.chdir('C:\\Users\\Rebecca\\Documents\\NU_MSPA\\thesis') # laptop

# load normalized FPKM table
FPKM_norm = pd.read_csv('data/fpkm_table_normalized.csv', 
                        index_col = 'gene_id \ rnaseq_profile_id')

# convert to numeric
FPKM_norm = FPKM_norm.astype(float)

# add column with mean expression level
FPKM_norm['mean_level'] = FPKM_norm.mean(axis = 1)

# reset index on dataframe
FPKM_norm = FPKM_norm.reset_index()

# nonzero mean indices
nonzero_mean = FPKM_norm.mean_level.nonzero()[0]

# Create dataframe for clustering: delete zero-mean rows and mean_level column
FPKM_clust = FPKM_norm.iloc[nonzero_mean].drop(['mean_level'], axis = 1)

# reset index to 'gene_id \ rnaseq_profile_id'
FPKM_clust = FPKM_clust.set_index('gene_id \ rnaseq_profile_id')

# heatmap
plt.figure(1)
heatmap_grid = sns.clustermap(FPKM_clust, 
                              method='complete', 
                              metric='correlation', 
                              figsize=(10,8),
                              cmap='PiYG')