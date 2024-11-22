# libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
})

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

filter_predictions_by_gene_peak <- function(method, starting_pred, cell_type, peaks, genes, out_file) {
	df_start = fread(starting_pred, sep="\t", header=TRUE)
	df_start$PeakName = with(df_start,  paste0(chr, "-", start, "-", end))
	df_filt = dplyr::filter(df_start, TargetGene %in% genes, PeakName %in% peaks)
	df_filt$CellType = cell_type

	pct_filtered = round(nrow(df_filt)/nrow(df_start) * 100, 2)
	message("Method: ", method, ", % rows retained: ", pct_filtered)

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



### MAIN
## files + params
# BMMC5
# dataset = "BMMC5_peaks_only"
# peaks = "/oak/stanford/groups/engreitz/Projects/scE2G/data/10x_BMMC_neurips2021_5_supergroups/detected_gene_peak/peaks_BMMC5"   # cell_type.bed (chr,start,end)
# genes = "/oak/stanford/groups/engreitz/Projects/scE2G/data/10x_BMMC_neurips2021_5_supergroups/detected_gene_peak/df.rna_detected_percent.5.bi.csv" # cols = gene, celltypes; rows = genes; boolean
# pred_dir ="/oak/stanford/groups/engreitz/Projects/scE2G/data/Genome_wide_1Mb_prediction/BMMC" # method/gene_peak_filtered/celltype.pairs.E2G.res.tsv.gz
# pred_subdir = "5_super_groups"
# cell_types = c("B", "Dendritic", "Erythroid", "Myeloid", "T") 

dataset = "BMMC22_peaks_only"
peaks = "/oak/stanford/groups/engreitz/Projects/scE2G/data/10x_BMMC_neurips2021/detected_gene_peak/peaks_BMMC22"   # cell_type.bed (chr,start,end)
genes = "/oak/stanford/groups/engreitz/Projects/scE2G/data/10x_BMMC_neurips2021/detected_gene_peak/df.rna_detected_percent.22.bi.csv" # cols = gene, celltypes; rows = genes; boolean
pred_dir ="/oak/stanford/groups/engreitz/Projects/scE2G/data/Genome_wide_1Mb_prediction/BMMC" # method/gene_peak_filtered/celltype.pairs.E2G.res.tsv.gz
pred_subdir = "22_cell_types"
cell_types = c("B1_B", "CD4_pos_T_activated", "CD4_pos_T_naive", "CD8_pos_T_naive", "CD8_pos_T", "CD14_pos_Mono", "CD16_pos_Mono",
	"cDC2", "Erythroblast", "G_M_prog", "HSC", "ID2_hi_myeloid_prog", "ILC", "Lymph_prog", "MK_E_prog", "Naive_CD20_pos_B", "NK", "Normoblast",
		"pDC", "Plasma_cell", "Proerythroblast", "Transitional_B") 


working_dir = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0623_CTS_predictions"
#methods = c("SnapATAC", "Signac", "Cicero", "FigR", "ScenicPlus", "DIRECTNET", "ArchR")
methods = c("Cicero")

## run
# read in expression matrices
genes = fread(genes, sep=",", header=TRUE)
colnames(genes)[1] = "TargetGene"

# iterate through cell types
for (i in 1:length(cell_types)){
	group = cell_types[i]
	this_genes = get_expressed_genes(genes, group) 
	this_peaks = get_accessible_peaks(peaks, group) 
	message("Cell type: ", group, ", # genes: ", length(this_genes), ", # peaks:  ", length(this_peaks))

	## iterate through methods
	for (m in 1:length(methods)){
		method = methods[m]
		if (method=="ScenicPlus") {file_name = "pairs.E2G.res.1000000.tsv.gz"} else {file_name = "pairs.E2G.res.tsv.gz"}
		starting_pred = file.path(pred_dir, method, pred_subdir, file_name)
		out_dir =  file.path(working_dir, dataset, method)
		if (i==1){ dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)}
		out_file = file.path(out_dir, paste0(group, ".pairs.E2G.res.tsv.gz"))

		#filter_predictions_by_gene_peak(method, starting_pred, group, this_peaks, this_genes, out_file)
		filter_predictions_by_peak_only(method, starting_pred, group, this_peaks, out_file)
	}
}