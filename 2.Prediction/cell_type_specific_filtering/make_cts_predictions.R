# libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
})

# helper function to parse GTF
extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- unlist(strsplit(gtf_attributes, " "))
  if(att_of_interest %in% att){
    return(gsub("\"|;","", att[which(att %in% att_of_interest)+1]))
  } else {
    return(NA)}
}

# map gene names from RNA matrix and expression metrics to gene reference used by scE2G via Ensembl ID
map_gene_names <- function(df_genes, gene_gtf_path, abc_genes_path){
	gene_ref <- fread(gene_gtf_path, header = FALSE, sep = "\t") %>%
		setNames(c("chr","source","type","start","end","score","strand","phase","attributes")) %>%
		dplyr::filter(type == "gene")
	gene_ref$gene_ref_name <- unlist(lapply(gene_ref$attributes, extract_attributes, "gene_name"))
	gene_ref$Ensembl_ID <- unlist(lapply(gene_ref$attributes, extract_attributes, "gene_id"))
	gene_ref <- dplyr::select(gene_ref, gene_ref_name, Ensembl_ID) %>%
		mutate(Ensembl_ID = sub("\\.\\d+$", "", Ensembl_ID)) %>% # remove decimal digits 
		distinct()
	
	abc_genes <- fread(abc_genes_path, col.names = c("chr", "start", "end", "name", "score", "strand", "Ensembl_ID", "gene_type")) %>%
		dplyr::select(name, Ensembl_ID) %>%
		rename(abc_name = name) %>%
		left_join(gene_ref, by = "Ensembl_ID") %>%
    	filter(!is.na(gene_ref_name)) %>%
		group_by(Ensembl_ID) %>% # remove cases where multiple genes map to one ensembl ID
		filter(n() == 1) %>%
		ungroup()

	gene_key <- abc_genes$abc_name
	names(gene_key) <- abc_genes$gene_ref_name # gencode_name = abc_name

	# remove genes not in our gene universe	
	df_genes <- df_genes %>% filter(TargetGene %in% names(gene_key)) %>% 
		mutate(TargetGene = gene_key[TargetGene])

	return(df_genes)
}

# return vector of expressed gene names
get_expressed_genes <- function(df_detected, cell_type) {
	df_genes = dplyr::select(df_detected, "TargetGene", starts_with(cell_type) & ends_with(cell_type))
	colnames(df_genes) = c('TargetGene', "expressed")
	df_genes = dplyr::filter(df_genes, expressed==TRUE)

	return(df_genes$TargetGene)
} 

# return vector of expressed peaks names (format: chr1-100004103-10000430)
get_accessible_peaks <- function(peak_dir, cell_type) {
	peak_file = file.path(peak_dir, paste0(cell_type, ".bed"))
	df_peaks = fread(peak_file, sep="\t", header=FALSE)
	colnames(df_peaks) = c("chr", "start", "end")

	df_peaks$PeakName = with(df_peaks, paste0(chr, "-", start, "-", end))

	return(df_peaks$PeakName)
}

# implement suggested threshold -> threshold_binary: Correlation >= 0.45 & FDR <= 1e-04 & VarQATAC > 0.25 & VarQRNA > 0.25
threshold_ArchR <- function(df, cell_type) {
	df$threshold_binary = ifelse(df$Correlation>=0.45 & df$FDR<=1e-04 & df$VarQATAC>0.25 & df$VarQRNA>0.25, 1, 0)

	pct_pos = round(sum(df$threshold_binary)/nrow(df)*100, 2)
	message("Biosample: ", cell_type, ", % positive: ", pct_pos)

	return(df)
}

# implement suggested threshold -> HC_binary: function_type==HC
threshold_DIRECTNET <- function(df, cell_type) {
	df$HC_binary = ifelse(df$function_type=="HC", 1, 0)

	pct_pos = round(sum(df$HC_binary)/nrow(df)*100, 2)
	message("Biosample: ", cell_type, ", % positive: ", pct_pos)

	return(df)
}

filter_predictions_by_gene_peak <- function(method, starting_pred, cell_type, peaks, genes, out_file) {
	df_start = fread(starting_pred, sep="\t", header=TRUE)
	df_start$PeakName = with(df_start,  paste0(chr, "-", start, "-", end))
	df_filt = dplyr::filter(df_start, TargetGene %in% genes, PeakName %in% peaks)
	df_filt$CellType = cell_type

	pct_filtered = round(nrow(df_filt)/nrow(df_start) * 100, 2)
	message("Method: ", method, ", % rows retained: ", pct_filtered)

	if (method == "DIRECTNET") {
		df_filt <- threshold_DIRECTNET(df_filt, cell_type)
	}
	 
	if (method == "ArchR") {
		df_filt <- threshold_ArchR(df_filt, cell_type)
	}

	fwrite(df_filt, out_file, col.names=TRUE, row.names=FALSE, sep="\t")
}

filter_predictions_by_peak_only <- function(method, starting_pred, cell_type, peaks, out_file) {
	df_start = fread(starting_pred, sep="\t", header=TRUE)
	df_start$PeakName = with(df_start,  paste0(chr, "-", start, "-", end))
	df_filt = dplyr::filter(df_start, PeakName %in% peaks)
	df_filt$CellType = cell_type

	pct_filtered = round(nrow(df_filt)/nrow(df_start) * 100, 2)
	message("Method: ", method, ", % rows retained: ", pct_filtered)

	fwrite(df_filt, out_file, col.names=TRUE, row.names=FALSE, sep="\t")
}

threshold_predictions_only <- function(method, starting_pred, cell_type, out_file) {
	df = fread(starting_pred, sep="\t", header=TRUE)
	df$CellType = cell_type

	if (method == "DIRECTNET") {
		df_thresh <- threshold_DIRECTNET(df, cell_type)
	}
	 
	if (method == "ArchR") {
		df_thresh <- threshold_ArchR(df, cell_type)
	}

	fwrite(df_thresh, out_file, col.names=TRUE, row.names=FALSE, sep="\t")

}





### MAIN
filter_to_cell_types <- function() {
	## files + params
	# dataset = "BMMC22"
	# peaks = "/oak/stanford/groups/engreitz/Projects/scE2G/data/10x_BMMC_neurips2021/detected_gene_peak/peaks_BMMC22"   # cell_type.bed (chr,start,end)
	# genes = "/oak/stanford/groups/engreitz/Projects/scE2G/data/10x_BMMC_neurips2021/detected_gene_peak/df.rna_detected_percent.22.bi.csv" # cols = gene, celltypes; rows = genes; boolean
	# cell_types = c("B1_B", "CD4_pos_T_activated", "CD4_pos_T_naive", "CD8_pos_T_naive", "CD8_pos_T", "CD14_pos_Mono", "CD16_pos_Mono",
	# 	"cDC2", "Erythroblast", "G_M_prog", "HSC", "ID2_hi_myeloid_prog", "ILC", "Lymph_prog", "MK_E_prog", "Naive_CD20_pos_B", "NK", "Normoblast",
	# 		"pDC", "Plasma_cell", "Proerythroblast", "Transitional_B") 

	# BMMC5
	# dataset = "BMMC5"
	# peaks = "/oak/stanford/groups/engreitz/Projects/scE2G/data/10x_BMMC_neurips2021_5_supergroups/detected_gene_peak/peaks_BMMC5"   # cell_type.bed (chr,start,end)
	# genes = "/oak/stanford/groups/engreitz/Projects/scE2G/data/10x_BMMC_neurips2021_5_supergroups/detected_gene_peak/df.rna_detected_percent.5.bi.csv" # cols = gene, celltypes; rows = genes; boolean
	# cell_types = c("B", "Dendritic", "Erythroid", "Myeloid", "T") 

	# PBMC9
	# dataset = "PBMC9"
	# peaks = "/oak/stanford/groups/engreitz/Projects/scE2G/data/Genome_wide_1Mb_prediction/PBMC9/detected_gene_peak/peaks_PBMC"   # cell_type.bed (chr,start,end)
	# genes = "/oak/stanford/groups/engreitz/Projects/scE2G/data/Genome_wide_1Mb_prediction/PBMC9/detected_gene_peak/df.rna_detected_percent.9.bi.csv" # cols = gene, celltypes; rows = genes; boolean
	# cell_types = c("10_cDC", "11_CD14.Mono.1", "12_CD14.Mono.2", "13_CD16.Mono", "17_B", "20_CD4.N1", "22_CD4.M", "24_CD8.CM", "25_NK") 

	# PBMC5
	# dataset = "PBMC5"
	# peaks = "/oak/stanford/groups/engreitz/Projects/scE2G/data/Genome_wide_1Mb_prediction/PBMC5/detected_gene_peak/peaks_PBMC"   # cell_type.bed (chr,start,end)
	# genes = "/oak/stanford/groups/engreitz/Projects/scE2G/data/Genome_wide_1Mb_prediction/PBMC5/detected_gene_peak/df.rna_detected_percent.5.bi.csv" # cols = gene, celltypes; rows = genes; boolean
	# cell_types = c("B", "DC", "Mono", "NK", "T") 

	gene_gtf_path <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_temp/scE2G/resources/genome_annotations/gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz"
	abc_genes_path <-"/oak/stanford/groups/engreitz/Users/sheth/scE2G_temp/scE2G/resources/genome_annotations/CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed"

	working_dir = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0623_CTS_predictions"
	pred_config_file = file.path(working_dir, "config", "CTS_pred_config.tsv")
	pred_config <- fread(pred_config_file) %>% filter(biosample == dataset)
	pred_key <- setNames(pred_config$pred_file, pred_config$method)

	methods = c("SnapATAC", "Signac", "Cicero", "FigR", "ScenicPlus", "DIRECTNET", "ArchR")

	## run
	# read in expression matrices
	genes = fread(genes, sep=",", header=TRUE)
	message("Genes before mapping names: ", nrow(genes))
	colnames(genes)[1] = "TargetGene"
	genes <- map_gene_names(genes, gene_gtf_path, abc_genes_path)
	message("Genes after mapping names: ", nrow(genes))

	# iterate through cell types
	for (i in 1:length(cell_types)){
		group = cell_types[i]
		this_genes = get_expressed_genes(genes, group)
		this_peaks = get_accessible_peaks(peaks, group) 
		message("Cell type: ", group, ", # genes: ", length(this_genes), ", # peaks:  ", length(this_peaks))

		## iterate through methods
		for (m in 1:length(methods)){
			method = methods[m]
			starting_pred = pred_key[method]

			out_dir = file.path(working_dir, dataset, method); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
			out_file = file.path(out_dir, paste0(group, ".pairs.E2G.res.tsv.gz"))

			if (method != "Cicero") {
				filter_predictions_by_gene_peak(method, starting_pred, group, this_peaks, this_genes, out_file)
			} else {
				filter_predictions_by_peak_only(method, starting_pred, group, this_peaks, out_file)
			}
		}
	}

}

just_apply_threshold <- function() {
	## just threshold K562 and GM12878 for DIRECTNET and ArchR
	working_dir = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0623_CTS_predictions"
	pred_config_file = file.path(working_dir, "config", "CTS_pred_config.tsv")
	pred_config <- fread(pred_config_file) %>% filter(biosample %in% c("K562_Xu", "GM12878"))
	print(pred_config)
	out_dir <- file.path(working_dir, "K562_GM12878"); dir.create(out_dir)

	for (i in 1:nrow(pred_config)) {
		print(i)
		starting_pred <- pred_config$pred_file[i]
		this_method <- pred_config$method[i]
		cell_type <- pred_config$biosample[i]

		out_dir_this <- file.path(out_dir, this_method); dir.create(out_dir_this)
		out_file_this <- file.path(out_dir_this, paste0(cell_type, ".pairs.E2G.res.tsv.gz"))

		threshold_predictions_only(this_method, starting_pred, cell_type, out_file_this)
	}
}

just_apply_threshold()
