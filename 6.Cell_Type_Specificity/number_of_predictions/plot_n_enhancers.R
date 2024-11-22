suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(tidyverse)
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(egg)
})

sort_samples <- function(cs_string){	
	cs <- strsplit(cs_string, ",") %>%  unlist() %>% sort() %>% paste(collapse = " | ")
	return(cs)

}

process_input_data <- function(input_file, input_id, out_dir) {
	res <- fread(input_file)
	print(sort_samples(res$split_samples[1]))
	res$split_samples <- lapply(res$split_samples, sort_samples) 
	res$scramble_samples <- lapply(res$scramble_samples, sort_samples)
	#	mutate(split_samples =sort_samples(split_samples),
	#		scramble_samples = sort_samples(split_samples)) %>%
	res <- res %>%	distinct() %>%
		group_by(supergroup, split_samples, scramble_samples) %>%
		summarize(n_total_split = mean(n_total_split), n_total_scramble = mean(n_total_scramble), n_shared = mean(n_shared)) %>%
		ungroup() %>% group_by(supergroup) %>%
		summarize(n_reps = n(), n_combos_split = n_distinct(split_samples), n_combos_scramble = n_distinct(scramble_samples),
			mean_total_split = mean(n_total_split), mean_total_scramble = mean(n_total_scramble), mean_shared = mean(n_shared),
			sd_total_split = sd(n_total_split), sd_total_scramble = sd(n_total_scramble),
			n_split_wins = sum(n_total_split > n_total_scramble)) %>%
		mutate(CI_high_total_scramble = 1.96 * sd_total_scramble + mean_total_scramble,
			CI_low_total_scramble = mean_total_scramble - 1.96 * sd_total_scramble,
			CI_high_total_split = 1.96 * sd_total_split + mean_total_split,
			CI_low_total_split = mean_total_split - 1.96 * sd_total_split)

	fwrite(res, file.path(out_dir, paste0("summary_", input_id, ".tsv")), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

	return(res)
} 

plot_n_enhancers <- function(res, clusters_per_sg, input_id, out_dir) {
	res <- dplyr::filter(res, supergroup != "BMMC5_Dendritic")
	res <- mutate(res, n_clusters = clusters_per_sg[supergroup],
		supergroup = paste0(gsub("BMMC5_", "", supergroup), "\n", n_clusters, " clusters"),
		fraction_wins = paste0(n_split_wins, "/", n_reps),
		max_CI = pmax(CI_high_total_split, CI_high_total_scramble)) 

	res_split <- dplyr::select(res, supergroup, mean_total = mean_total_split, CI_low_total = CI_low_total_split, CI_high_total = CI_high_total_split) %>%
		mutate(Group = "Granular clusters", color = "#792374")
	res_scramble <- dplyr::select(res, supergroup, mean_total = mean_total_scramble, CI_low_total = CI_low_total_scramble, CI_high_total = CI_high_total_scramble) %>%
		mutate(Group = "Supergroups", color = "#96a0b3")
	res_all <- rbind(res_split, res_scramble)

	res_text <- group_by(res, supergroup, fraction_wins) %>% 
		summarize(max_CI = max(max_CI)) %>% 
		mutate(Group = "None") %>%
		arrange(max_CI)
	cp <- c(`Granular clusters`="#792374", `Supergroups`="#96a0b3")

	res_all$supergroup <- factor(res_all$supergroup, levels = res_text$supergroup, ordered = TRUE)

	bar <- ggplot(res_all, aes(y = mean_total, x = supergroup, fill = Group)) +
		geom_col(position = position_dodge(0.9), width = 0.9) +
		geom_errorbar(aes(ymin=CI_low_total, ymax=CI_high_total),
			position = position_dodge(0.9), width=0.2, linewidth = 0.25, color = "#000000") +
		geom_text(data = res_text, aes(x = supergroup, y = (max_CI + 5000), label = fraction_wins), color = "#000000", size = 2.5) +
		scale_fill_manual(values = cp) +
		xlab("Cell type supergroup") + ylab("Number of unique enhancer-gene links") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position = "right")

	ggsave(file.path(out_dir, paste0("bar_plot_v2_", input_id, ".pdf")), bar, height = 3.5, width = 4.75)
}


input_id = "100r_17s_3d"
input_file = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0916_global_properties/n_enhancers_BMMC/n_enhancers_100r_17s_3d.tsv"
# supergroup	rep	split_samples	scramble_samples	n_total_split	n_total_scramble	n_shared	n_uniq_split	n_uniq_scramble
out_dir <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0916_global_properties/n_enhancers_BMMC"
sample_key_file <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0916_global_properties/workflow/n_enhancer_analysis/BMMC_split_scramble_sample_key.tsv"
clusters_per_sg <- c(BMMC5_B=3, BMMC5_T=4, BMMC5_Myeloid=3, BMMC5_Dendritic=2, BMMC5_Erythroid=3)

res <- process_input_data(input_file, input_id, out_dir)
plot_n_enhancers(res, clusters_per_sg, input_id, out_dir)
