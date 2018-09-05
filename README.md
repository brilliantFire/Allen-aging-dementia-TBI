# Allen-aging-dementia-TBI
Repo for loading and analyzing the data from the Allen Institute of Brain Science's [Aging, Dementia, &amp; TBI Study](http://aging.brain-map.org/). Inputs and outputs of scripts in `ETL` and `EDA` are stored in `data`. The `work` folder contains things currently being worked on.  

The `blog-notebooks` folder has some Jupyter notebooks I made in the process of blogging about this project. See [my Jekyll blog](http://blog.vislaywade.com/).  

This project is being used to fulfill the requirements of the Master of Science in Predictive Analytics program from Northwestern University.  Thanks for checking it out!

#### ../ETL
* *'expected_count_TPM_FPKM_load.py/.R'* - Loads data from study website (R and Python versions)  
* *'subset_count_data.R'* - Subsets RNA-seq counts by brain region or donor sex   

#### ../EDA
* *'plotly_mds_plots_html.R'* - Generates [HTMLwidget multidimensional scaling plots](http://blog.vislaywade.com/interactive-MDS-plots-w-plotly/) of RNA-seq data using the `plotly` library  
* *'differential_expression_dispersions_exact_tests.R'* - hypothesis testing for differential gene expression  
* *'differential_expression_visualizations.R'* - visualizations & analysis of differential analysis results
