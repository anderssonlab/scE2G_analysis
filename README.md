## scE2G manuscript analyses

Scripts and notebooks for data preprocessing and analyses performed in Sheth, Qiu, et al. 2024 (add link)

* For the scE2G tool, please see to the scE2G model and pipeline [`repo`](https://github.com/EngreitzLab/sc-E2G)
* For the scE2G prediction results, please see https://data.broadinstitute.org/scE2Gpredictions/ 
* For the benchmarking pipelines, please see to the individual benchmarking repos for [`CRISPR`](https://github.com/EngreitzLab/CRISPR_comparison/tree/main), [`eQTL`](https://github.com/EngreitzLab/eQTLEnrichment/tree/integrated), and [`GWAS`](https://github.com/EngreitzLab/GWAS_E2G_benchmarking)

## Contents:

### 1.Preprocess_Data
Code for generating input files for E2G prediction models

### 2.Prediction
Code for applying published models for E2G prediction 

### 3.Benchmarking
Code and results for CRISPR, eQTL and GWAS benchmarking analyses
* eQTL
	* `2024_1004_hd_pred_fine_eQTL`: results for "fine-grained" eQTL analysis
	* `2024_1003_scE2G_eQTLCat_supergroups`: results for "coarse-grained" eQTL analysis
	* `eQTL_Catalogue_v7`: code to collate and group fine-mapped eQTLs from the eQTL Catalogue
		* `config/dataset_metadata_subset.tsv`: Subset of metadata from https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/data_tables/dataset_metadata.tsv with condition_label == "naive" and quant_method == "ge". The columns `tissue_fine` and `tissue_coarse` contain assigments of datasets used for the "fine-grained" and "coarse-grained" eQTL analyses, respectively. All other columns are from the source.
* GWAS
	* `2024_1006_hd_clusters`: results for GWAS benchmarking of enhancer predictions in fine-grained cell types
	* `2024_1005_supergroups`: results for GWAS benchmarking of enhancer predictions in cell type supergroups, K562, and GM12878

### 4.Robustness_Evaluation
Code for evaluation of model robustness basing on downsampling data

### 5. Predictions_Properties
Code for evaluation of properties of scE2G predictions

### 6.Cell_Type_Specificity
Code for evaluation of cell-type specific E2G predictions

### 7.GWAS_Variant_Interpretation
Code for analyses related to gene prioritization and interpreation of GWAS loci 
