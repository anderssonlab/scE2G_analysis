suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(tidyverse)
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(cowplot)
  library(ggpubr)
  library(ggdist)
})


get_distinct_combinations_sorted <- function(df, col1, col2) {
	 df <- df %>%
		rename(s1 = !!sym(col1), s2 = !!sym(col2)) %>% # Temporarily rename for processing
		mutate(
		s1 = as.character(s1), s2 = as.character(s2), 
		s1_sorted = pmin(s1, s2), s2_sorted = pmax(s1, s2) # Sort columns
		) %>%
		filter(s1_sorted != s2_sorted) %>% # Exclude rows where s1 == s2
		select(-s1, -s2) %>% # Drop temporary columns
		rename(!!sym(col1) := s1_sorted, !!sym(col2) := s2_sorted) %>% # Restore original column names
		distinct(!!sym(col1), !!sym(col2), .keep_all = TRUE) %>% # Keep unique combinations
		mutate(pair_id = paste0(!!sym(col1), "_", !!sym(col2))) # Add pair ID
  
  return(df)
}

process_all_cluster_input <- function(input_file, sample_key_file){
	sample_key <- fread(sample_key_file) # biosample	supergroup	cell_type	dataset
	sample_key_A <- setNames(sample_key, paste0(colnames(sample_key), "A"))
	sample_key_B <- setNames(sample_key, paste0(colnames(sample_key), "B"))

	res_uniq <- fread(input_file) %>%
		dplyr::filter(biosampleA != biosampleB) %>%
		rename(biosample1 = biosampleA, biosample2 = biosampleB) %>%
		mutate(biosampleA = pmin(biosample1, biosample2), biosampleB = pmax(biosample1, biosample2),
			jaccard = nSharedPredAwB / (nTotalPredA + nTotalPredB - nSharedPredAwB),
			pair_id = paste0(biosampleA, "_", biosampleB)) %>%
		select(-c(biosample1, biosample2, nSharedPredAwB, nTotalPredA, nTotalPredB)) %>%
  		distinct(biosampleA, biosampleB, .keep_all = TRUE) %>%
		dplyr::filter(biosampleA %in% sample_key$biosample, biosampleB %in% sample_key$biosample) %>%
		left_join(sample_key_A, by = "biosampleA") %>% left_join(sample_key_B, by = "biosampleB") 
	#print(res_uniq)

	# within dataset group
	win1 <- expand.grid(biosampleA = sample_key$biosample[sample_key$dataset == "BMMC"], biosampleB = sample_key$biosample[sample_key$dataset == "BMMC"])
	win2 <- expand.grid(biosampleA = sample_key$biosample[sample_key$dataset == "PBMC"], biosampleB = sample_key$biosample[sample_key$dataset == "PBMC"])
	win_pairs<- rbind(win1, win2) %>% get_distinct_combinations_sorted("biosampleA", "biosampleB")
	win <- dplyr::filter(res_uniq, pair_id %in% win_pairs$pair_id) %>%
		mutate(group = "Same dataset")
	#print(win)

	# between datasets
	btw_pairs <- expand.grid(biosampleA = sample_key$biosample[sample_key$dataset == "BMMC"],
		biosampleB = sample_key$biosample[sample_key$dataset == "PBMC"]) %>%
		get_distinct_combinations_sorted("biosampleA", "biosampleB")
	btw <- dplyr::filter(res_uniq, pair_id %in% btw_pairs$pair_id) %>%
		mutate(group = "Different dataset")
	#print(btw)

	# sg 
	sg <- dplyr::filter(res_uniq, datasetA == datasetB, supergroupA == supergroupB) %>%
		mutate(group = "Same cell type, same dataset")

	# ct
	ct <- dplyr::filter(res_uniq, datasetA != datasetB, cell_typeA==cell_typeB) %>%
		mutate(group = "Same cell type, different dataset")

	full <- rbind(win, btw, sg, ct)
}

summarize_groups <- function(res, out_dir) {
	res_summ <- group_by(res, group) %>%
		summarize(n_pairs = n(), mean_jaccard = mean(jaccard), mean_Pearson = mean(Pearson),
			sd_jaccard = sd(jaccard), sd_Pearson = sd(Pearson)) %>%
		mutate(CI_jaccard_low = mean_jaccard - 1.96*sd_jaccard, CI_jaccard_high = mean_jaccard + 1.96*sd_jaccard,
			CI_Pearson_low = mean_Pearson - 1.96*sd_Pearson, CI_Pearson_high = mean_Pearson + 1.96*sd_Pearson)
}

pairwise_comparisons <- function(res, metric, out_dir) {
	combos = expand.grid(group1 = unique(res$group), group2 = unique(res$group)) %>%
		dplyr::filter(group1 != group2) %>%
		get_distinct_combinations_sorted("group1", "group2")
	print(combos)

	# iterate through rows
	res_list <- vector("list", nrow(combos))
  	for (i in 1:nrow(combos)) {
		group1 = combos$group1[i]
		group2 = combos$group2[i]

		group1_val <- dplyr::filter(res, group == group1) %>% pull(!!sym(metric))
		group2_val <- dplyr::filter(res, group == group2) %>% pull(!!sym(metric))

		# message(group1, ": ", length(group1_val), " values")
		# message(group2, ": ", length(group2_val), " values")

		ttest <- t.test(x = group1_val, y = group2_val, alternative = "two.sided", paired = FALSE)
		 
		 res_list[[i]] <- as.data.frame(list(group1 = group1, group2 = group2, metric = metric,
		 	group1_mean = mean(group1_val), group2_mean = mean(group2_val),
			delta_group1_minus_groups2 = mean(group1_val) - mean(group2_val), p = ttest$p.value))
  	}
	
	res <- rbindlist(res_list) %>% as_tibble() %>%
		mutate(p_adjust = p.adjust(p, method = "bonferroni"),
			significant = p_adjust < 0.05) %>%
		distinct() %>%
		arrange(p_adjust)

	fwrite(res, file.path(out_dir, paste0(metric, "_significance_table.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
	return(res)

}

plot_cluster_grouping_comparisons <- function(cluster_res, j_signif, p_signif, out_dir) {
	cp_j <- c(`Same cell type, same dataset`="#00488d", `Same dataset`="#5496ce", `Same cell type, different dataset`="#006eae", `Different dataset`="#9bcae9")
	cp_p <- c(`Same cell type, same dataset`="#9b241c", `Same dataset`="#dc6464", `Same cell type, different dataset`="#c5373d", `Different dataset`="#e9a0a5")
	summ <- summarize_groups(cluster_res, out_dir)
	
	n_pairs_key <- setNames(summ$n_pairs, summ$group)
	label_key <- c(`Same cell type, same dataset`="Same cell type group,\nsame dataset", `Same dataset`="Same dataset",
		`Same cell type, different dataset`="Same cell type group,\ndifferent dataset", `Different dataset`="Different dataset")

	cluster_res <- mutate(cluster_res,
			n_pairs = n_pairs_key[group], jaccard_col = cp_j[group], Pearson_col = cp_p[group],
			n_pairs_text = paste0(n_pairs, " pairs"),
			group_label = paste0(label_key[group], "\n", n_pairs_text),
			group = factor(group, levels = names(cp_j), ordered = TRUE)) %>%
		arrange(group) 
	cp_j <- setNames(unique(cluster_res$jaccard_col), unique(cluster_res$group_label))
	cp_p <- setNames(unique(cluster_res$Pearson_col), unique(cluster_res$group_label))
	cluster_res$group_label = factor(cluster_res$group_label, levels = names(cp_j), ordered = TRUE)

	j <- ggplot(cluster_res, aes(x = group_label, y = jaccard)) +
		stat_eye(aes(fill = group_label), side = "both", shape = 16, point_size = 2, slab_linewidth = 0) +
		scale_fill_manual(values = cp_j) +
		geom_signif(xmin = c(1,3), xmax = c(2,4),
		 	annotations = c("***", "***"), y_position = c(0.58, 0.58), tip_length = 0.02,
			size = 0.25, extend_line = -0.02, vjust = 0.7, textsize = 2.5) +
		scale_y_continuous(limits = c(0, 0.7), breaks = c(0, 0.2, 0.4, 0.6)) +
		labs(x = "Category of cell type pair", y = "Jaccard similarity") +
		theme_classic() + theme(axis.text = element_text(size = 6, color = "#000000"), axis.title = element_text(size = 8),
			axis.text.x = element_blank(), axis.title.x = element_blank(),
			axis.ticks = element_line(color = "#000000"), legend.position = "none")

	p <- ggplot(cluster_res, aes(x = group_label, y =Pearson)) +
		stat_eye(aes(fill = group_label), side = "both", shape = 16, point_size = 2, slab_linewidth = 0) +
		scale_fill_manual(values = cp_p) +
		geom_signif(xmin = c(1,3), xmax = c(2,4),
		 	annotations = c("***", "***"), y_position = c(0.67, 0.67), tip_length = 0.02,
			size = 0.25, extend_line = -0.02, vjust = 0.7, textsize = 2.5) +
		scale_y_continuous(limits = c(0, 0.7), breaks = c(0, 0.2, 0.4, 0.6)) +
		labs(x = "Category of cell type pair", y = "Pearson correlation") +
		theme_classic() + theme(axis.text = element_text(size = 6, color = "#000000"), axis.title = element_text(size = 8),
			axis.ticks = element_line(color = "#000000"), legend.position = "none")

	
	grid <- plot_grid(j, p, nrow = 2, axis = "lr", align = "hv")
	ggsave2(file.path(out_dir, paste0("jaccard_pearson_between_groups.pdf")), grid, height = 4.5, width = 3.5)
	# fwrite(res, file.path(out_dir, paste0(corr_col, "_summarized_across_BMMC.tsv")), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}

cluster_res_file <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0214_new_global_properties/enhancer_similarity/correlation_across_clusters.tsv"
cluster_sample_key_file <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0214_new_global_properties/config/cell_type_groups.tsv"

out_dir <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0214_new_global_properties/enhancer_similarity"
correlation_cols <- c("jaccard", "Pearson", "Pearson_log1p", "Spearman") #"Pearson_log1p" #"Pearson_log1p" #\tPearson\tPearson_log1p""

dir.create(out_dir, showWarnings = FALSE)
res <- process_all_cluster_input(cluster_res_file, cluster_sample_key_file)
j_signif <- pairwise_comparisons(res, "jaccard", out_dir)
p_signif <- pairwise_comparisons(res, "Pearson", out_dir)
plot_cluster_grouping_comparisons(res, j_signif, p_signif, out_dir)
