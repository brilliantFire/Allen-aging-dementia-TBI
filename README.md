# Allen-aging-dementia-TBI
Repo for loading and analyzing the data from the Allen Institute of Brain Science's [Aging, Dementia, &amp; TBI Study](http://aging.brain-map.org/). The goal of the project is to develop features from gene expression levels for over 50,000 genes to be used along with neuropathological measurements, demographic information, and medical history data from 107 donors in models of dementia probability. Most scripts were written using [Anaconda](https://www.anaconda.com/) distributions of R 3.5.1 and [Jupyter Lab](http://jupyter.org/). There are a few Python versions that I've played with but the bulk of the project is in R.  

Inputs and outputs of scripts in `ETL`, `EDA` , `feature-engineering`, `final-dataset-prep`, and `lasso-models-and-evaluation` are stored in `data`. The `work` folder contains scripts in progress. 

The `blog-notebooks` folder has some Jupyter notebooks I made in the process of blogging about this project. See [my Jekyll blog](http://blog.vislaywade.com/).  

This project is being used to fulfill the requirements of the Master of Science in Predictive Analytics program (now [Master's in Data Science](https://sps.northwestern.edu/masters/data-science/)) from Northwestern University.  Thanks for checking it out!

#### ../ETL
* *`expected_count_TPM_FPKM_load.py/.R`* - Loads data from study website (R and Python versions)  
* *`subset_count_data.R`* - Subsets RNA-seq counts by brain region or donor sex   

#### ../EDA
* *`plotly_mds_plots_html.R`* - Generates [HTMLwidget multidimensional scaling plots](http://blog.vislaywade.com/interactive-MDS-plots-w-plotly/) of RNA-seq data using the `plotly` library  
* *`differential_expression_dispersions_exact_tests.R`* - hypothesis testing for differential gene expression  
* *`differential_expression_visualizations.R`* - visualizations & analysis of differential analysis results

#### ../feature-engineering
* *`clustering_experiments_gene_expression.R`* - Uses `clValid` to test different clustering methods and parameters on gene expression levels data
* *`feature_engineering_dataset_construction`* - Runs k-medoids clustering via [CLARA](http://www.sthda.com/english/articles/27-partitioning-clustering-essentials/89-clara-clustering-large-applications/); adds medoids as features to tables of neuropathological, demographic, and medical history data for each patient/donor; accesses the [mygene.info](http://mygene.info/) database for gene ontologies and creates cloud visualizations of biological process terms for each cluster

#### ../final-dataset-prep  
* *`data_preparation.R`* - Preparation of final dataset for modeling, including: dropping variables, converting categorical features to factor variables & reordering levels, and imputation of missing variables by CART (classification & regression trees) using the `mice` library; also includes location metric (mean & median) comparisons for all numeric variables

#### ../lasso-models-and-evaluation  
* *`lasso_models_and_evaluation.R`* - [LASSO](https://en.wikipedia.org/wiki/Lasso_(statistics)) (Least Absolute Shrinkage and Selection Operator)-penalized logistic regression using the `glmnet` library; Labeling at two different cutoffs probability; Model evaluation using confusion matrices and [ROC](https://en.wikipedia.org/wiki/Receiver_operating_characteristic) curves

##### This README last updated on: 2018-10-20
