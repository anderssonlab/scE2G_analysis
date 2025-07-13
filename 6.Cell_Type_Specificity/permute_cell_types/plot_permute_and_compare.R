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

plot_n_enhancers <- function(res_summ, res_sign, out_dir) {
	label_key <- c(L3_real_total = "True cell type\nlabels", L3_from_L2_total = "Permuted labels\nwithin supergroup", L3_from_L1_total = "Permuted labels\nacross dataset")
	cp <- c(`True cell type\nlabels` = "#430b4e", `Permuted labels\nwithin supergroup` = "#a64791", `Permuted labels\nacross dataset` = "#96a0b3")

	res <- res_summ %>% mutate(sample_label = label_key[sample_type]) %>% 
		arrange(mean_enhancers) %>% 
		mutate(sample_label = factor(sample_label, levels = label_key, ordered = TRUE),
			mean_enhancers = mean_enhancers/1e3, CI_low = CI_low/1e3, CI_high = CI_high/1e3)

	res_sign <- res_sign %>% 
		mutate(group1 = label_key[group1], group2 = label_key[group2],
			y_loc = c(925, 950, 925))
	
	print(res_sign)
	
	## graphing params
	y_max <- 1e6 / 1e3
	x_lab <- paste0("Sample category\n(N = ", mean(res$n_reps), " iterations)")

	p <- ggplot(res, aes(y = mean_enhancers, x = sample_label)) +
		geom_col(aes(fill = sample_label)) +
		geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width=0.2, linewidth = 0.25, color = "#000000") +
		geom_signif(xmin = res_sign$group1, xmax = res_sign$group2,
			annotations = res_sign$p_stars, y_position = res_sign$y_loc, tip_length = 100,
			size = 0.25, extend_line = -10, vjust = 0.7, textsize = 2.5) +
		scale_fill_manual(values = cp) +
		ylim(c(0, y_max)) +
		labs(x = x_lab, y = "# unique enhancer-gene links across dataset (x 10^3)") +
		coord_flip() +
		theme_classic() + theme(axis.text = element_text(size = 7, color = "#000000"),
			axis.title = element_text(size = 8),
			axis.ticks = element_line(color = "#000000"), legend.position = "none")

	ggsave(file.path(out_dir, "L1_totals.pdf"), p, height = 1.6, width = 4)
}

pairwise_comparisons <- function(res_all, res_summ) {
	res_all <- select(res_all, -iter)
	groups = data.frame(t(combn(c(colnames(res_all)), 2)))
  	colnames(groups) = c("group1", "group2")

  # iterate through rows
	res_list <- vector("list", nrow(groups))
  	for (i in 1:nrow(groups)) {
		group1 = groups$group1[i]
		group2 = groups$group2[i]

		group1_val <- res_all %>% pull(!!sym(group1))
		group2_val <- res_all %>% pull(!!sym(group2))

		ttest <- t.test(x = group1_val, y = group2_val, alternative = "greater", paired = TRUE)
		 
		res_list[[i]] <- as.data.frame(list(group1 = group1, group2 = group2,
		 	group1_mean = mean(group1_val), group2_mean = mean(group2_val),
			delta_group1_minus_groups2 = mean(group1_val) - mean(group2_val),
			p = ttest$p.value))
  	}
	
	res <- rbindlist(res_list) %>% as_tibble() %>%
		mutate(p_stars = case_when(p < 0.001 ~ "***",
									p < 0.01 ~ "**",
									p < 0.05 ~ "*",
									TRUE ~ ""))
	return(res)
}

process_inputs <- function(input_files, out_dir) {
	# iter	L3_real_total	L3_from_L2_total	L3_from_L1_total

	res_all <- lapply(input_files, fread) %>% 
		rbindlist() %>% as.data.frame()

	res_summ <- res_all %>%
		pivot_longer(cols = -iter, names_to = "sample_type", values_to = "total_links") %>%
		group_by(sample_type) %>%
		summarize(n_reps = n(),
			mean_enhancers = mean(total_links),
			sd_enhancers = sd(total_links)) %>%
		mutate(CI_high = mean_enhancers + 1.96 * sd_enhancers,
			CI_low = mean_enhancers - 1.96 * sd_enhancers)
	fwrite(res_summ, file.path(out_dir, "L1_total_summary.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

	sign_table <- pairwise_comparisons(res_all, res_summ)
	fwrite(sign_table, file.path(out_dir, "significance_table.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

	return(list(summary = res_summ, signif = sign_table))
}

res_dir <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_permute_cell_types/scE2G/results"
res_ids <- c("2025_0220_iter10_seed13", "2025_0220_iter10_seed19")
out_dir <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0214_new_global_properties/permute_cell_types"; dir.create(out_dir)

input_files <- file.path(res_dir, res_ids, "L1_totals.tsv")

res_list <- process_inputs(input_files, out_dir)
res_summ <- res_list[["summary"]]
res_sign <- res_list[["signif"]]
plot_n_enhancers(res_summ, res_sign, out_dir)
