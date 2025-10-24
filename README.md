## scE2G manuscript analyses

Scripts, notebooks, and results for data preprocessing and analyses performed in Sheth, Qiu, _et al._ 2024, [`preprint`](https://www.biorxiv.org/content/10.1101/2024.11.23.624931v1)

* For the scE2G tool, please see to the scE2G model and pipeline [`repo`](https://github.com/EngreitzLab/scE2G)
* For the scE2G prediction results, please see https://data.broadinstitute.org/scE2Gpredictions/ 
* For the benchmarking pipelines, please see to the individual benchmarking repos for [`CRISPR`](https://github.com/EngreitzLab/CRISPR_comparison/tree/main), [`eQTL`](https://github.com/EngreitzLab/eQTLEnrichment/tree/integrated), and [`GWAS`](https://github.com/EngreitzLab/GWAS_E2G_benchmarking)

<hr>

## System requirements

### Hardware requirements

For the analyses, a standard computer equipped with 32+ GB RAM is necessary.

### Software requirements

The analyses were conducted on Linux systems.

The versions of the software tested are listed in the Methods section of Sheth, Qiu, _et al._ 2024

<hr>

## Contents:

### 1.Preprocess_Data
Code for generating input files for E2G prediction models. Each directory corresponds to a dataset listed in Table S10

### 2.Prediction
Code for applying published models for E2G prediction. The file structure is organized as [dataset]/[method]

### 3.Benchmarking
Code and results for CRISPR, eQTL and GWAS benchmarking analyses
* CRISPR
	* `Fig.2/K562_Xu/sc.250425_2.K562.Xu.cv`: benchmarking with training CRISPR data < 1 Mb for Xu _et al._ K562 dataset (Fig. 2b,c)
	* `Fig.2/K562_Xu/sc.250425.K562.Xu.cv.ignore_TPM`: benchmarking with training CRISPR data < 1 Mb for Xu _et al._ K562 dataset, including scE2G<sup>multiome</sup> predictions that ignore gene expression filtering (Fig. S6a)
	* `Fig.2/K562_Wang/sc.250425.K562.Wang.ignore_TPM`: benchmarking with training CRISPR data < 1 Mb for Wang _et al._ K562 dataset, including scE2G<sup>multiome</sup> predictions that ignore gene expression filtering (Fig. S6b)
	* `Fig.2/K562_This_Study/sc.250425.K562.IGVF.ignore_TPM`: benchmarking with training CRISPR data < 1 Mb for K562 dataset in this study, including scE2G<sup>multiome</sup> predictions that ignore gene expression filtering (Fig. S6c)
	* `Fig.2/K562_Xu_all_distances/sc.250425_2.K562.Xu.cv.5M`: benchmarking with training CRISPR data covering all distances for Xu _et al._ K562 dataset (Fig. S6d)
	* `Fig.2/K562_Xu_default/sc.250425.K562.default.Xu.cv`: benchmarking with training CRISPR data < 1 Mb for Xu _et al._ K562 dataset, using existing models with their default settings (Fig. S7)
	* `Fig.3/sc.250425.K562.SEACell90_PBMC`: benchmarking with training CRISPR data < 1 Mb for Xu _et al._ K562 dataset, using three approaches to calculating Kendall correlation and corresponding scE2G models (Fig. 3c)
	* `heldout_crispr_benchmark`: benchmarking with held-out CRISPR data < 1 Mb for K562, GM12878, HCT116, Jurkat, and WTC11 datasets (Fig. 2d)
* eQTL
	* `2024_1004_hd_pred_fine_eQTL`: results for "fine-grained" eQTL analysis
	* `2024_1003_scE2G_eQTLCat_supergroups`: results for "coarse-grained" eQTL analysis
	* `eQTL_Catalogue_v7`: code to collate and group fine-mapped eQTLs from the eQTL Catalogue
		* `config/dataset_metadata_subset.tsv`: Subset of metadata from https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/data_tables/dataset_metadata.tsv with condition_label == "naive" and quant_method == "ge". The columns `tissue_fine` and `tissue_coarse` contain assigments of datasets used for the "fine-grained" and "coarse-grained" eQTL analyses, respectively. All other columns are from the source.
* GWAS
	* `2024_1006_hd_clusters`: results for GWAS benchmarking of enhancer predictions in fine-grained cell types
	* `2024_1005_supergroups`: results for GWAS benchmarking of enhancer predictions in cell type supergroups, K562, and GM12878

### 4.Model_Interpretation
Code and results from feature analysis of scE2G, and code for evaluation of model robustness basing on downsampling data

### 5. Predictions_Properties
Code for evaluation of properties of scE2G predictions

### 6.Cell_Type_Specificity
Code and results for evaluation of cell-type specific E2G predictions

### 7.GWAS_Variant_Interpretation
Code for analyses related to gene prioritization and interpreation of GWAS loci 
