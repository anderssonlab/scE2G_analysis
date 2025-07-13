suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(ggdist)
  library(cowplot)
})

get_distribution_within_cell_type <- function(n, key, metric, model_params, include_promoters) {
	cluster <- key$biosample[n]
	score_threshold <- model_params["score_threshold"] %>% as.numeric()
	model_name <- model_params["model_name"]
	score_column <- model_params["score_column"]
	pred_file <- key[[model_name]][n]

	message(model_name, ", ", cluster, " (", n, "/", nrow(key), ")")

	# filter predictions
	pred <- fread(pred_file)
	pred <- dplyr::filter(pred, !!sym(score_column) >= score_threshold)
	if (!include_promoters) {
		pred <- pred %>% dplyr::filter(class != "promoter")
	}

	# get values
	if (metric == "enh_region_size") {
		res <- dplyr::select(pred, chr, start, end) %>% distinct() %>%
			mutate(value = end - start)

	} else if (metric == "log10_dist_to_tss") {
		res <- dplyr::select(pred, chr, start, end, TargetGene, distance) %>% distinct() %>%
			mutate(value = log10(distance))

	} else if (metric == "dist_to_tss") {
		res <- dplyr::select(pred, chr, start, end, TargetGene, distance) %>% distinct() %>%
			mutate(value = distance)

	} else if (metric == "num_enh_per_gene") {
		res <- dplyr::select(pred, chr, start, end, TargetGene) %>% distinct() %>%
			group_by(TargetGene) %>%
			tally() %>%
			mutate(value = n)

	} else if (metric == "num_gene_per_enh") {
		res <- dplyr::select(pred, chr, start, end, TargetGene) %>% distinct() %>%
			group_by(chr, start, end) %>%
			tally() %>%
			mutate(value = n)

	} else {
		message("Please chose valid metric.")
	}

	# add other info
	res <- dplyr::select(res, value) %>%
		mutate(model = model_name, cluster = cluster)

	return(res)
}

plot_distribution_one_model <- function(res, other_res, model_name, metric) {
	# plotting params
	n_elements <- paste0("# elements: ", nrow(res))
	n_clusters <- paste0("# clusters: ", length(unique(res$cluster)))
	mean <- paste0("Mean: ", round(mean(res$value), 2))
	median <- paste0("median: ", round(median(res$value), 2))
	min <- paste0("min: ", round(min(res$value), 2))
	max <- paste0("max: ", round(max(res$value), 2))
	stats_label <- paste0(n_clusters, ", ", n_elements, "\n", mean, ", ", median, ", ", min, ", ", max)
	message(stats_label)
	color <- ifelse(model_name == "multiome_powerlaw_v3", "#792374", "#006479")

	discrete_metrics <-  c("num_enh_per_gene", "num_gene_per_enh", "enh_region_size")

	combined_vals <- c(res$value, other_res$value)
	quant_high <- quantile(combined_vals, 0.9)[1]; message(quant_high)
	quant_low <- quantile(combined_vals, 0.1)[1]; message(quant_low)
	
	# optionally, set axis limits and force data to this range
	if (metric == "dist_to_tss") {
		axis_lim <- c(750, 250000)
	} else if (metric == "enh_region_size") {
		axis_lim <- c(500, 1200)
	} else if (metric == "num_enh_per_gene") {
		axis_lim <- c(0, 20)
	} else if (metric == "num_gene_per_enh") {
		axis_lim <- c(0, 6)
	} else if (metric %in% c("none")) {  #right skew (set lower lim)
		axis_lim <- c(round(quant_low, 2), max(combined_vals))
	} else if (metric %in% c("dist_to_tss", "log10_dist_to_tss", "enh_region_size", "num_enh_per_gene")) { # left skew (set upper lim)
		axis_lim <- c(min(combined_vals), round(quant_high, 2))
	} else {
		axis_lim <- c(min(combined_vals), max(combined_vals))
	}

	axis_lim_label <- paste0("Axis limits (for reference): ", axis_lim[1], "-", axis_lim[2])
	message(axis_lim_label)
	res$value <- pmin(res$value, axis_lim[2])
	res$value <- pmax(res$value, axis_lim[1])

	label_key <- c(enh_region_size = "Width of enhancer element (bp)",
		log10_dist_to_tss = "log10(Distance to target gene TSS (bp))",
		dist_to_tss = "Distance to target gene TSS (bp)",
		num_enh_per_gene = "Number of enhancers per gene",
		num_gene_per_enh = "Number of genes per enhancer")
	
	if (metric %in% discrete_metrics) {
		if (metric == "enh_region_size") {
			bin_width =50
		} else {
			bin_width = 1
		}
		x <- ggplot(res, aes(x = value)) +
				stat_slabinterval(side = "top", shape = 16, point_size = 2, slab_linewidth = 0,
					slab_fill = color, alpha = 0.9, density = "histogram",
					breaks = breaks_fixed(width = bin_width), align = align_boundary(at = 0)) + 
				labs(x = label_key[metric], y = NULL, title = stats_label, caption = axis_lim_label) +
				xlim(axis_lim) +
				theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
					plot.title = element_text(size = 8), legend.position = "None", plot.caption = element_text(size = 7))
	} else {
		x <- ggplot(res, aes(x = value)) +
				stat_slabinterval(side = "top", shape = 16, point_size = 2, slab_linewidth = 0,
					slab_fill = color, alpha = 0.9, density = "histogram") + 
				labs(x = label_key[metric], y = NULL, title = stats_label, caption = axis_lim_label) +
				xlim(axis_lim) +
				theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
					plot.title = element_text(size = 8), legend.position = "None", plot.caption = element_text(size = 7))
	}

	if (metric %in% c("dist_to_tss")) { # scale 
		x <- x + scale_x_log10()
	}

	return(x)
}

plot_distribution_both_models <- function(metric, multiome_model_params, atac_model_params, ds_include, ct_exclude, include_promoters, out_dir) {
	if (metric == "dist_to_tss") {
		rds_multiome_file <- file.path(out_dir, paste0("log10_", metric, "_multiome.rds"))
		rds_atac_file <- file.path(out_dir, paste0("log10_", metric, "_atac.rds"))
	} else {
		rds_multiome_file <- file.path(out_dir, paste0(metric, "_multiome.rds"))
		rds_atac_file <- file.path(out_dir, paste0(metric, "_atac.rds"))
	}

	if (file.exists(rds_multiome_file)) {
		multiome_res <- readRDS(rds_multiome_file)
	} else {
		multiome_key <- fread(multiome_model_params["sample_key_file"]) %>%
			dplyr::filter(dataset %in% ds_include, !(biosample %in% ct_exclude))
		multiome_res <- lapply(1:nrow(multiome_key), get_distribution_within_cell_type, key = multiome_key, metric = metric,
			model_params = multiome_model_params, include_promoters = include_promoters) %>%
			rbindlist() %>% as.data.frame()
		saveRDS(multiome_res, rds_multiome_file)
	}
	
	if (file.exists(rds_atac_file)) {
		atac_res <- readRDS(rds_atac_file)
	} else {
		atac_key <- fread(atac_model_params["sample_key_file"]) %>%
			dplyr::filter(dataset %in% ds_include, !(biosample %in% ct_exclude))
		atac_res <- lapply(1:nrow(atac_key), get_distribution_within_cell_type, key = atac_key, metric = metric,
			model_params = atac_model_params, include_promoters = include_promoters) %>%
			rbindlist() %>% as.data.frame()
		saveRDS(atac_res, rds_atac_file)
	}

	if (metric == "dist_to_tss") {
		multiome_res <- multiome_res %>% mutate(value = 10 ** value)
		atac_res <- atac_res %>% mutate(value = 10 ** value)
	}
	m_plot <- plot_distribution_one_model(multiome_res, atac_res, multiome_model_params["model_name"], metric)
	a_plot <- plot_distribution_one_model(atac_res, multiome_res, atac_model_params["model_name"], metric)

	# assemble
	gr <- plot_grid(m_plot, a_plot, nrow = 2, ncol = 1, rel_heights = c(1, 1), align = "hv", axis = "lr")
	return(gr)

}

##### RUN
sample_key <- fread("/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0214_new_global_properties/config/pred_sample_key.tsv")
#sample_key <- fread("/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0916_global_properties/config/atac_sample_key.tsv")

if (TRUE){
	# params
	out_dir <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0214_new_global_properties/prediction_properties/distributions"
	dir.create(out_dir, showWarnings = FALSE)
	multiome_model_params <- c(model_name = "multiome_powerlaw_v3", score_threshold = 0.177, score_column = "E2G.Score.qnorm",
		sample_key_file = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0214_new_global_properties/config/pred_sample_key.tsv" )
	atac_model_params <- c(model_name = "scATAC_powerlaw_v3", score_threshold = 0.174, score_column = "E2G.Score.qnorm",
		sample_key_file = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0214_new_global_properties/config/pred_sample_key.tsv")
	include_promoters <- FALSE
	ds_include <- c("BMMC22", "PBMC9", "Islets")
	#ds_include <- c("GM12878")
	ct_exclude <- c("BMMC22_ID2_hi_myeloid_prog")



	# plot metrics across all granular cell types
	#metrics <- c("enh_region_size",  "log10_dist_to_tss", "num_enh_per_gene", "num_gene_per_enh")
	#metrics <- c("num_gene_per_enh", "num_enh_per_gene", "enh_region_size", "dist_to_tss")
	metrics <- c("dist_to_tss")
	for (m in metrics) {
		p <- plot_distribution_both_models(m, multiome_model_params, atac_model_params, ds_include, ct_exclude, include_promoters, out_dir)
		ggsave2(file.path(out_dir, paste0(m, "_granular_clusters_2M_noKG.pdf")), p, width=4, height=4)
		# ggsave2(file.path(out_dir, paste0(m, "_GM12878_test.pdf")), p, width=4, height=4)
	}

}

