suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(ggdist)
  library(ggpubr)
  library(cowplot)
})

make_enr_recall_table <- function(result_dir, key_file, category_group, biosamples_exclude) {
	enr_file <- file.path("/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-integrated/eQTLEnrichment/results", result_dir, "scE2G_multiome", "enrichmentTables", "enrichmentTable.0to30000Kb.tsv")
	recall_file <- file.path("/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-integrated/eQTLEnrichment/results", result_dir, "scE2G_multiome", "recallTables", "recallTable.byDistance.tsv")

	# enr file: Biosample	GTExTissue	nVariantsOverlappingEnhancers	nVariantsGTExTissue	nCommonVariantsOverlappingEnhancers	nCommonVariants	enrichment	CI_enr_low	CI_enr_high	SE_log_enr	p_adjust_enr	method	distance_min	distance_max
	# recall file (need to filter) GTExTissue	nVariantGenePairsOverlappingEnhancers	nVariantsOverlappingEnhancersCorrectGene	Biosample	distance_min	distance_max	total.variants	recall.total	recall.linking	correctGene.ifOverlap	method
	enr <- fread(enr_file)
	recall <- fread(recall_file)
	key <- fread(key_file)
	clust_map <- dplyr::select(key, Biosample, clustering) %>% distinct()

	if (category_group == "any_exact") {
		category_values <- c("exact", "almost_exact")
	} else if (category_group == "not_unmatched") {
		category_values <-  c("exact", "almost_exact", "similar")
	}

	res <- left_join(enr, recall, by=c("GTExTissue", "Biosample", "distance_min", "distance_max", "method")) %>%
		dplyr::filter(nVariantsGTExTissue >= 35,!is.na(enrichment), !is.na(recall.linking)) %>%
		dplyr::filter(!(Biosample %in% biosamples_exclude)) %>%
		left_join(key, by = c("GTExTissue", "Biosample")) %>%
		dplyr::select(-clustering) %>%
		mutate(category = replace_na(category, "unmatched"),
			category = ifelse(category %in%  category_values, category_group, category),
			category = factor(category, levels = c("exact", "almost_exact", "any_exact", "similar", "not_unmatched", "unmatched"), ordered = TRUE)) %>%
		left_join(clust_map, by = "Biosample") %>% 
		dplyr::filter(!is.na(clustering))

	#test <- dplyr::filter(res, is.na(clustering)) %>% dplyr::select(Biosample, GTExTissue, clustering, category)
	#print(test)

	return(res)
}

plot_enrichment_recall_scatter <- function(er, out_dir) {
	no_facet <- dplyr::select(er, enrichment, recall.linking)
	group_counts <- group_by(er, clustering, category) %>%
		summarize(n = n()) %>%
		mutate(count_label = paste0("N = ", n))

	# facet by clustering and category
	s <- ggplot(er, aes(x = recall.linking, y = enrichment)) +
		facet_grid(rows = vars(clustering), cols = vars(category), scales = "fixed", axes = "all") +
		geom_point(data = no_facet, shape = 16, size = 1, alpha = 0.5, color = "#c5cad7") +
		geom_point(shape = 16, size = 2, alpha = 0.5, color = "#1c2a43") +
		geom_text(data = group_counts, aes(x = 0.2, y = 2, label = count_label)) +
		scale_y_log10() + 
		labs(x = "Recall (Fraction of eQTLs in enhancer linked to eGene)", y = "log10 Enrichment (eQTLs versus common variants)") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
			strip.background = element_blank(), aspect.ratio=1)

	ggsave(file.path(out_dir, "enrichment_recall_scatter_v0.pdf"), s, width = 8, height = 6)
		
}

plot_enrichment_recall_distributions <- function(er, out_dir) {
	group_counts <- group_by(er, clustering, category) %>%
		summarize(n = n(),
		max_enr = max(enrichment),
		max_rec_linking = max(recall.linking),
		avg_enr = mean(enrichment),
		avg_rec_linking = mean(recall.linking)) %>%
		mutate(count_label = paste0("N = ", n),
			grouping = paste0(clustering, "\n", category)) %>%
		arrange(avg_enr)

	er$grouping = paste0(er$clustering, "\n", er$category)
	er$grouping = factor(er$grouping, levels = group_counts$grouping, ordered = TRUE)

	# facet by clustering and category
	e <- ggplot(er, aes(x = enrichment, y = grouping)) +
		geom_vline(xintercept = 1, linetype = "dashed", color = "#c5cad7") +
		stat_eye(side = "both", shape = 16, point_size = 2, slab_linewidth = 0, slab_fill = "#96a0b3") +
		geom_text(data = group_counts, aes(x = max_enr + 3, y = grouping, label = count_label), size = 2.5) +
		scale_x_log10() + 
		labs(y = "Prediction cluster\nBiosample matching", x = "log10 Enrichment (eQTLs versus common variants)") +
		theme_classic() + theme( axis.text = element_text(size = 7), axis.title = element_text(size = 8))
	ggsave(file.path(out_dir, "enrichment_distributions.pdf"), e, width = 6, height = 4)

	group_counts <- arrange(group_counts, avg_rec_linking)
	er$grouping = factor(er$grouping, levels = group_counts$grouping, ordered = TRUE)
	r <- ggplot(er, aes(x = recall.linking, y = grouping)) +
		stat_eye(side = "both", shape = 16, point_size = 2, slab_linewidth = 0, slab_fill = "#c5cad7") +
		geom_text(data = group_counts, aes(x = max_rec_linking + 0.03, y = grouping, label = count_label), size = 2.5) +
		labs(y = "Prediction cluster\nBiosample matching", x = "Recall (Fraction of eQTLs in enhancer linked to eGene)") +
		theme_classic() + theme( axis.text = element_text(size = 7), axis.title = element_text(size = 8))
	ggsave(file.path(out_dir, "recall_distributions.pdf"), r, width = 6, height = 4)
}

compare_metric_across_groups <- function(er, metric, out_dir) {
	er$grouping = paste0(er$clustering, "_", er$category)
	
	groups = data.frame(t(combn(c(unique(er$grouping)), 2)))
  	colnames(groups) = c("group1", "group2")

  # iterate through rows
	res_list <- vector("list", nrow(groups))
  	for (i in 1:nrow(groups)) {
		group1 = groups$group1[i]
		group2 = groups$group2[i]

		group1_val <- dplyr::filter(er, grouping == group1) %>% pull(!!sym(metric))
		group2_val <- dplyr::filter(er, grouping == group2) %>% pull(!!sym(metric))

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

plot_main_comparison <- function(er, enr_sign, rec_sign, out_dir) {
	groups <- c("granular_any_exact", "granular_similar", "granular_not_unmatched", "granular_unmatched")
	
	# format for plotting
	er_groups <- group_by(er, category) %>%
		summarize(n_points = n(), mean_enr = mean(enrichment), mean_rec = mean(recall.linking)) %>%
		mutate(category_name = case_when(category == "any_exact" ~ "Same cell type",
																		category == "not_unmatched" ~ "Same cell type\nor lineage",
																		category == "similar" ~ "Same lineage",
																		category == "unmatched" ~ "Unmatched"),
					category_alpha = case_when(category == "any_exact" ~ 0.75,
																		category == "not_unmatched" ~ 0.75,
																		category == "similar" ~ 0.75,
																		category == "unmatched" ~ 0.33),
					category_color = case_when(category == "any_exact" ~ "#430b4e",
																		category == "not_unmatched" ~ "#430b4e",
																		category == "similar" ~ "#a64791",
																		category == "unmatched" ~ "#96a0b3"),
					category_label = paste0(category_name, "\n(N = ", n_points, ")"),
					clustering = "granular", group = paste0(clustering, "_", category)) %>%
		arrange(mean_enr)

	cat_key <- setNames(er_groups$category_label, er_groups$group)
	enr_rank <- setNames(rank(er_groups$mean_enr), er_groups$group[rank(er_groups$mean_enr)])
	rec_rank <- setNames(rank(er_groups$mean_rec), er_groups$group[rank(er_groups$mean_rec)])
	cp <- setNames(er_groups$category_color, er_groups$category_label)
	
	enr_lim <- c(0, 15.5)

	enr_sign <- dplyr::filter(enr_sign, group1 %in% groups, group2 %in% groups, significant == TRUE) %>%
		mutate(stars = case_when(p_adjust < 0.001 ~ "***", p_adjust < 0.01 ~ "**", p_adjust < 0.05 ~ "*", TRUE ~ ""),
			category1 = cat_key[group1], category2 = cat_key[group2],
			loc1 = enr_rank[group1], loc2 = enr_rank[group2],
			y_position = ifelse(loc1+loc2==4, enr_lim[2]-0.5, enr_lim[2]-1))
	print(enr_sign %>% select(loc1, loc2, y_position))

	rec_sign <- dplyr::filter(rec_sign, group1 %in% groups, group2 %in% groups, significant == TRUE) %>%
		mutate(stars = case_when(p_adjust < 0.001 ~ "***", p_adjust < 0.01 ~ "**", p_adjust < 0.05 ~ "*", TRUE ~ ""),
			category1 = cat_key[group1], category2 = cat_key[group2],
			loc1 = rec_rank[group1], loc2 = rec_rank[group2],
			y_position = ifelse(loc1+loc2==4, max(er$recall.linking)+ 0.018, max(er$recall.linking) + 0.01))

	er <- left_join(er, er_groups, by = "category") %>%
		mutate(category_label = factor(category_label, levels = (cat_key[enr_rank]), ordered = TRUE)) %>%
		arrange(category_label)

	s <- ggplot(er, aes(x = recall.linking, y = enrichment, color = category_label, alpha = category_alpha)) +
		geom_point(shape = 16, size = 2) +
		ylim(enr_lim) +
		scale_color_manual(values = cp, name = "Prediction cell type/\neQTL biosample category") +
		scale_alpha_identity() +
		#scale_y_log10() + 
		labs(x = "Recall\nFraction of eQTLs in enhancer linked to correct eGene", y = "Enrichment\neQTLs versus common variants") +
		theme_classic() + theme( axis.text = element_text(size = 7), axis.title = element_text(size = 8), axis.ticks = element_line(color = "#000000"), legend.position = "right",
			legend.text = element_text(size = 7), legend.title = element_text(size = 8), aspect.ratio=1)
	ggsave(file.path(out_dir, "granular_scatter.pdf"), s, height = 4, width = 4.5)

	er$category_label <- factor(er$category_label, levels = cat_key[names(enr_rank)], ordered = TRUE)
	e <- ggplot(er, aes(x = category_label, y = enrichment)) +
		stat_eye(aes(fill = category_label), side = "both", shape = 16, point_size = 2, slab_linewidth = 0) +
		ylim(enr_lim) +
		scale_fill_manual(values = cp, name = "Prediction cell type/eQTL biosample category") +
		geom_signif(xmin = enr_sign$loc1, xmax = enr_sign$loc2,
			annotations = enr_sign$stars, y_position = enr_sign$y_position, tip_length = 0.01,
			size = 0.25, extend_line = -0.02, vjust = 0.7, textsize = 2.5) +
		labs(x = "Prediction cluster/eQTL biosample category", y = "Enrichment\neQTLs versus common variants") +
		#coord_flip() +
		theme_classic() + theme(axis.text = element_text(size = 7, color = "#000000"), axis.title = element_text(size = 8),
			axis.ticks = element_line(color = "#000000"), legend.position = "none")

	er$category_label <- factor(er$category_label, levels = cat_key[names(rec_rank)], ordered = TRUE)
	r <- ggplot(er, aes(x = category_label, y = recall.linking)) +
		stat_eye(aes(fill = category_label), side = "both", shape = 16, point_size = 2, slab_linewidth = 0) +
		#ylim(enr_lim) +
		scale_fill_manual(values = cp, name = "Prediction cell type/eQTL biosample category") +
		geom_signif(xmin = rec_sign$loc1, xmax = rec_sign$loc2,
			annotations = rec_sign$stars, y_position = rec_sign$y_position, tip_length = 0.01,
			size = 0.25, extend_line = -0.02, vjust = 0.7, textsize = 2.5) +
		labs(x = "Prediction cluster/eQTL biosample category", y = "Recall\nFraction of eQTLs in enhancer linked to correct eGene") +
		#coord_flip() +
		theme_classic() + theme( axis.text = element_text(size = 7, color = "#000000"), axis.title = element_text(size = 8),
			axis.ticks = element_line(color = "#000000"), legend.position = "none")

		
	grid <- cowplot::plot_grid(e,r, nrow = 1, ncol = 2, align = "hv", rel_widths = c(1, 1))
	ggsave2(file.path(out_dir, "granular_violins_v.pdf"), grid, height = 4, width = 6)

	grid2 <- cowplot::plot_grid(s + theme(legend.position = "none"), e, r, nrow = 1, ncol = 3, align = "vh", rel_widths = c(0.85, 0.7, 0.7))
	# ggsave2(file.path(out_dir, "granular_scatter_and_violins.pdf"), grid2, height = 3.5, width = 8.5)

}

## INPUT FILES
# sg_enr_file = "/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-integrated/eQTLEnrichment/results/2024_1109_sg_pred_fine_eQTL/scE2G_multiome/enrichmentTables/enrichmentTable.0to30000Kb.tsv"
# sg_recall_file = "/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-integrated/eQTLEnrichment/results/2024_1109_sg_pred_fine_eQTL/scE2G_multiome/recallTables/recallTable.byDistance.tsv"
# sg_key_file ="/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0623_CTS_predictions/config/sg_pred_fine_eQTL_key.tsv"

# hd_enr_file = "/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-integrated/eQTLEnrichment/results/2024_1004_hd_pred_fine_eQTL/scE2G_multiome/enrichmentTables/enrichmentTable.0to30000Kb.tsv"
# hd_recall_file = "/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-integrated/eQTLEnrichment/results/2024_1004_hd_pred_fine_eQTL/scE2G_multiome/recallTables/recallTable.byDistance.tsv"
# hd_key_file ="/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0623_CTS_predictions/config/hd_pred_fine_eQTL_key.tsv"

result_dirs <- c("2025_0224_scE2G_only")
key_file <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0623_CTS_predictions/config/biosample_matching_hd_pred_fine_eQTL.tsv"
biosamples_exclude <- c("Islets_endothelial", "Islets_immune", "Islets_mesenchymal")
out_dir = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0623_CTS_predictions/eQTL_enrichment_specificity"; dir.create(out_dir)

er <- lapply(result_dirs, make_enr_recall_table, key_file, "any_exact", biosamples_exclude) %>%
	rbindlist() %>% distinct(GTExTissue, Biosample, .keep_all = TRUE) %>% as_tibble()
fwrite(er, file.path(out_dir, "annotated_enrichment_recall.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

message("Unique prediction cell types: ", length(unique(er$Biosample)))
message("Unique eQTL biosamples: ", length(unique(er$GTExTissue)))


plot_enrichment_recall_scatter(er, out_dir)
plot_enrichment_recall_distributions(er, out_dir)

enr_sign <- compare_metric_across_groups(er, "enrichment", out_dir)
rec_sign <- compare_metric_across_groups(er, "recall.linking", out_dir)

plot_main_comparison(er %>% dplyr::filter(clustering == "granular"), enr_sign, rec_sign, out_dir)

