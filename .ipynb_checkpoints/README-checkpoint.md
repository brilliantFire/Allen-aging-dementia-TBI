# Allen-aging-dementia-TBI
Repo for loading and analyzing the data from the Allen Institute of Brain Science's [Aging, Dementia, &amp; TBI Study](http://aging.brain-map.org/). Inputs and outputs of scripts in `ETL`, `EDA` , and `feature-engineering` are stored in `data`. The `work` folder contains scripts in progress.  

The `blog-notebooks` folder has some Jupyter notebooks I made in the process of blogging about this project. See [my Jekyll blog](http://blog.vislaywade.com/).  

This project is being used to fulfill the requirements of the Master of Science in Predictive Analytics program from Northwestern University.  Thanks for checking it out!

#### ../ETL
* *'expected_count_TPM_FPKM_load.py/.R'* - Loads data from study website (R and Python versions)  
* *'subset_count_data.R'* - Subsets RNA-seq counts by brain region or donor sex   

#### ../EDA
* *'plotly_mds_plots_html.R'* - Generates [HTMLwidget multidimensional scaling plots](http://blog.vislaywade.com/interactive-MDS-plots-w-plotly/) of RNA-seq data using the `plotly` library  
* *'differential_expression_dispersions_exact_tests.R'* - hypothesis testing for differential gene expression  
* *'differential_expression_visualizations.R'* - visualizations & analysis of differential analysis results

#### ../feature-engineering
* *`clustering_experiments_gene_expression.R`* - Uses `clValid` to test different clustering methods and parameters on gene expression levels data
* *`feature_engineering_dataset_construction`* - Runs k-medoids clustering via [CLARA](http://www.sthda.com/english/articles/27-partitioning-clustering-essentials/89-clara-clustering-large-applications/); adds medoids as features to tables of neuropathological, demographic, and medical history data for each patient/donor; accesses the [mygene.info](http://mygene.info/) database for gene ontologies and creates cloud visualizations of biological process terms for each cluster

##### This README last updated on: 2018-09-18