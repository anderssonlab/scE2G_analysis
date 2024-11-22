suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
  library(optparse)
})

option_list <- list(
    make_option(c("--supergroup"), type="character"),
	make_option(c("--dataset"), type="character"),
    make_option(c("--seed"), type="numeric"),
    make_option(c("--delta"), type="numeric"))
  
opt = parse_args(OptionParser(option_list=option_list))

sample_key_file <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0916_global_properties/workflow/n_enhancer_analysis/BMMC_split_scramble_sample_key.tsv"
clusters_per_sg <- c(BMMC5_B=3, BMMC5_T=4, BMMC5_Myeloid=3, BMMC5_Dendritic=2, BMMC5_Erythroid=3)

key_df <- fread(sample_key_file) %>%
	dplyr::filter(dataset == opt$dataset, supergroup == opt$supergroup)
file_key <- setNames(key_df$reformatted_pred, key_df$sample)

if (opt$dataset == "BMMC_scramble") {
	set.seed(opt$seed * opt$delta)
	n_samples <- nrow(key)
	sample_names <- sample(names(file_key), size = clusters_per_sg[opt$supergroup], replace = FALSE)

} else if (opt$dataset == "BMMC_split") {
	set.seed(opt$seed)
	sample_names <- group_by(key_df, cluster) %>%
		summarize(selected = sample(sample, 1)) %>%
		pull(selected)
}

# output
file_paths <- file_key[sample_names]

# output file paths as a single line (space-separated)
cat(paste(file_paths, collapse = " "), "\n")

# output sample names as a single line (comma-separated)
cat(paste(sample_names, collapse = ","), "\n")