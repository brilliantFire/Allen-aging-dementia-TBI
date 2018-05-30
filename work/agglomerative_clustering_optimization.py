# -*- coding: utf-8 -*-
"""
~*~*~*~* Agglomerative Clustering Optimization *~*~*~*~
Rebecca Vislay Wade
20 Apr 2018
--------------------------------------------------------------------------------
Performs a massive optimization for hierarchical clustering on a reduced version
of the normalized FPKM table from the Aging, Dementia, and TBI Study 
(http://aging.brain-map.org/). The results are saved to .csv and plotted.

The loop constructs clustering solutions for n_clusters values from 20 to 500
in increments of 20.

*** THIS TAKES A LONG TIME TO RUN AS WRITTEN ***
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

from sklearn.neighbors import kneighbors_graph
connectivity_matrix = kneighbors_graph(FPKM_clust, 
                                       n_neighbors = 20, 
                                       mode = 'connectivity', 
                                       include_self = True, 
                                       n_jobs = -1)
###
# list of n_cluster values between 20 and 500 in increments of 20 to try
clusters = list(range(20, 500, 20))

# try all three linkage types
linkages = ['ward', 'complete', 'average']

# loop consructs models, fits them, calculates CH_score, and puts the score 
# into a DataFrame for each type of linkage
random.seed(420)

# Initialize dataframe
results_df = pd.DataFrame(clusters, columns =['N_Clusters'])
#results_df.head()

for linkage_type in linkages:
    # Initialize results list
    results = []
    for num_clusters in clusters:
        model = AgglomerativeClustering(n_clusters=num_clusters,
                                        connectivity = connectivity_matrix,
                                        linkage = linkage_type)
        model.fit(FPKM_clust)
        score = calinski_harabaz_score(FPKM_clust, model.labels_)
        # Append scores to results list
        results.append(score)       
    # add results list as new column in DataFrame
    results_df[linkage_type]=results

# Save results_df to file
results_df.to_csv('thesis_code/output/csv/agglomerative_clustering_optimization_results_2018-04-23.csv')

# 3-panel figure CH_Score versus N_Clusters for each linkage type
plt.figure().set_size_inches(10,8)
for linkage_type in linkages:
    plt.subplot(3,1,(linkages.index(linkage_type)+1))
    plt.plot('N_Clusters', linkage_type, data=results_df, marker='o', color='darkgreen')
    plt.title(linkage_type)
plt.tight_layout()
plt.show()















