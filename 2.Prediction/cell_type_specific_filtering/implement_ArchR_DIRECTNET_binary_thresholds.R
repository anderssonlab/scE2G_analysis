# libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
})

# implement suggested threshold -> threshold_binary: Correlation >= 0.45 & FDR <= 1e-04 & VarQATAC > 0.25 & VarQRNA > 0.25
threshold_ArchR <- function(biosample, starting_pred, out_file) {
	df = fread(starting_pred, sep="\t", header=TRUE)
	df$threshold_binary = ifelse(df$Correlation>=0.45 & df$FDR<=1e-04 & df$VarQATAC>0.25 & df$VarQRNA>0.25, 1, 0)

	pct_pos = round(sum(df$threshold_binary)/nrow(df)*100, 2)
	message("Biosample: ", biosample, ", % positive: ", pct_pos)

	fwrite(df, out_file, col.names=TRUE, row.names=FALSE, sep="\t")
}

# implement suggested threshold -> HC_binary: function_type==HC
threshold_DIRECTNET <- function(biosample, starting_pred, out_file) {
	df = fread(starting_pred, sep="\t", header=TRUE)
	df$HC_binary = ifelse(df$function_type=="HC", 1, 0)

	pct_pos = round(sum(df$HC_binary)/nrow(df)*100, 2)
	message("Biosample: ", biosample, ", % positive: ", pct_pos)

	fwrite(df, out_file, col.names=TRUE, row.names=FALSE, sep="\t")
}

### MAIN
dataset = "BMMC"
working_dir = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0623_CTS_predictions"
sample_key_file = file.path(working_dir, "thresholding_sample_key.tsv")
out_dir = file.path(working_dir, "binarized", dataset) # .../binarized/dataset/method/cell_type.pairs.E2G.res.tsv.gz

sample_key = fread(sample_key_file, header=TRUE)
sample_key$ArchR_new = ""
sample_key$DIRECTNET_new = ""

## filter ArchR
dir.create(file.path(out_dir, "ArchR"), recursive=TRUE, showWarnings=FALSE)
message("ArchR")
for (i in 1:nrow(sample_key)){
	biosample = sample_key$biosample[i]
	starting_pred = sample_key$ArchR[i]
	out_file = file.path(out_dir, "ArchR", paste0(biosample, ".pairs.E2G.res.tsv.gz"))
	sample_key$ArchR_new[i] = out_file

	threshold_ArchR(biosample, starting_pred, out_file)
}

dir.create(file.path(out_dir, "DIRECTNET"), recursive=TRUE, showWarnings=FALSE)
message("DIRECTNET")
for (i in 1:nrow(sample_key)){
	biosample = sample_key$biosample[i]
	starting_pred = sample_key$DIRECTNET[i]
	out_file = file.path(out_dir, "DIRECTNET", paste0(biosample, ".pairs.E2G.res.tsv.gz"))
	sample_key$DIRECTNET_new[i] = out_file

	threshold_DIRECTNET(biosample, starting_pred, out_file)
}

fwrite(sample_key, sample_key_file, row.names=FALSE, col.names=TRUE, sep="\t")




