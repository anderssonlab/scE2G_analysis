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

plot_scatter <- function(df, x_metric, y_metric, label_key, dataset_key) {
	# format data
	n_label <- paste0("N = ", nrow(df))
	x_label <- paste0(label_key[x_metric], "\n", n_label)
	vline_loc <- ifelse(x_metric == "log10_unique_fragments", log10(2e6), mean(df[[y_metric]]))
	vline_alpha <- ifelse(x_metric == "log10_unique_fragments", 1, 0)

	# model name colors
	model_cp <- c(multiome_powerlaw_v3 = "#792374", scATAC_powerlaw_v3 = "#006479")

	# scatter plot
	s <- ggplot(df, aes(x = !!sym(x_metric), y = !!sym(y_metric), alpha = cluster_alpha)) +
		geom_vline(xintercept = vline_loc, alpha = vline_alpha, linetype = "dashed", color = "#96a0b3") +
		geom_point(aes(color = model_name), shape = 16, size = 1.5) +
		labs(x = x_label, y = label_key[y_metric]) +
		scale_color_manual(values = model_cp) +
		scale_alpha_identity() +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), 
			axis.ticks = element_line(color = "#000000"), aspect.ratio = 1, legend.position = "None")

	# margin distributions
	# x <- ggplot(df, aes(y = 1, x = !!sym(x_metric))) +
	# 	stat_halfeye(side = "top", fill = "#96a0b3", color = "#000000") + 
	# 	labs(x = NULL, y = NULL) +
	# 	theme_classic() + theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())

	# y <- ggplot(df, aes(x = 1, y = !!sym(y_metric))) +
	# 	stat_halfeye(side = "right", fill = "#96a0b3", color = "#000000") + 
	# 	labs(x = NULL, y = NULL) +
	# 	theme_classic() + theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())

	# assemble
	# gr <- plot_grid(x, NULL, s, y, nrow = 2, ncol = 2, rel_widths = c(5, 2), rel_heights = c(2, 5), align = "hv", axis = "lrtb")
	# return(gr)
	return(s)
}

plot_distribution <- function(df, metric, label_key) {
	# plotting params
	n_label <- paste0("N = ", nrow(df))
	model_cp <- c(multiome_powerlaw_v3 = "#792374", scATAC_powerlaw_v3 = "#006479")
	x_label <- paste0(label_key[metric], "\n", n_label)
	
	x <- ggplot(df, aes(x = !!sym(metric), y = model_name, color = model_name, alpha = cluster_alpha)) +
		geom_dots(side = "top", layout = "swarm", linewidth = 0, shape = 16, size = 1.5) + 
		labs(x = x_label, y = NULL) +
		scale_color_manual(values = model_cp) +
		scale_alpha_identity() +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
			axis.ticks = element_line(color = "#000000"), 
			axis.text.y = element_blank(), legend.position = "None")
}

plot_scatter_simple <- function(df, x_metric, y_metric, label_key) {
	# format data
	n_label <- paste0("N = ", nrow(df))
	x_label <- paste0(label_key[x_metric], "\n", n_label)



	# scatter plot
	s <- ggplot(df, aes(x = !!sym(x_metric), y = !!sym(y_metric), alpha = cluster_alpha)) +
		geom_point(aes(color = model_name), shape = 16, size = 1.5, color = "#1c2a43") +
		labs(x = x_label, y = label_key[y_metric]) + scale_alpha_identity() +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), aspect.ratio = 1,
				axis.ticks = element_line(color = "#000000"), legend.position = "None")
}

##### MAIN

## get inputs from snakemake
# sample_key <- fread(snakemake@input$sample_key)
# score_threshold <- snakemake@params$score_threshold %>% as.numeric()
# score_column <- snakemake@params$score_column
# include_promoters <- ifelse(snakemake@params$include_promoters %in% c("TRUE", "True"), TRUE, FALSE) 
stats_file <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0214_new_global_properties/prediction_properties/scE2G_prediction_properties.tsv"
out_dir <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0214_new_global_properties/prediction_properties"
sample_key_file <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0214_new_global_properties/config/pred_sample_key.tsv"

# labels 
metric_names <- c("fragments_total", "cell_count", "umi_count", "n_enh_elements", "n_genes_with_enh", "n_enh_gene_links",
		"mean_genes_per_enh", "mean_enh_per_gene", "mean_dist_to_tss", "mean_log10_dist_to_tss", "mean_enh_width",
		"n_genes_active_promoter", "n_genes_not_expressed", "log10_unique_fragments", "n_frag_per_cell", "log10_cell_count", "log10_umi_count")
metric_labels <- c("# unique ATAC fragments", "# cells", "# RNA UMIs", "# enhancer elements", "# genes with 1+ enhancer", "# enhancer-gene links",
	"Mean # genes per enhancer", "Mean # enhancers per gene", "Mean distance to TSS (bp)", "Mean log10(distance to TSS) (bp)",  "Mean enhancer width (bp)",
	"# genes considered", "# genes not expressed", "log10(# unique ATAC fragments)", "Mean # ATAC fragments per cell", "log10(# cells)", "log10(# RNA UMIs)")
label_key <- setNames(metric_labels, metric_names) # name: label


# read in and format data
sample_key <- fread(sample_key_file) %>% 
	select(biosample, biosample_name, dataset)

stats <- fread(stats_file) %>% 
	left_join(sample_key, by=c("cluster" = "biosample"))

# filter to dataset
ds_include <- c("BMMC22", "PBMC9", "Islets")
ct_exclude <- c("BMMC22_ID2_hi_myeloid_prog")
stats_for_summ <- stats %>%
	dplyr::filter(dataset %in% ds_include, !(cluster %in% ct_exclude), model_name == "multiome_powerlaw_v3") 

calcs <- stats_for_summ %>% 
	mutate(total_enh_size = n_enh_elements * mean_enh_width,
		total_enh_log10_dist = n_enh_gene_links * log10(mean_dist_to_tss))
mean_size <- sum(calcs$total_enh_size/sum(calcs$n_enh_elements))
mean_dist <- sum(calcs$total_enh_log10_dist/sum(calcs$n_enh_gene_links))
message("mean enhancer size: ", mean_size)
message("total number of enhancers: ", sum(calcs$n_enh_elements))
message("mean enhancer log10 distance to tss: ", mean_dist)
message("total number of links: ", sum(calcs$n_enh_gene_links))

summary <- stats_for_summ %>% pivot_longer(cols = -c(cluster, model_name, dataset, biosample_name), names_to = "metric", values_to = "value") %>% 
	group_by(metric) %>%
	summarize(mean = mean(value), median = median(value), min = min(value), max = max(value), n_clusters = n())
fwrite(summary, file.path(out_dir, "multiome_stats_summary_granular_2M_noKG.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# pivot wider
df <- stats %>% 
	mutate(num_unique_fragments_M = fragments_total/1e6,
		log10_unique_fragments = log10(fragments_total),
		log10_cell_count = log10(cell_count),
		n_frag_per_cell = fragments_total/cell_count,
		mean_log10_dist_to_tss = log10(mean_dist_to_tss),
		cluster_alpha = ifelse((dataset %in% ds_include & !(cluster %in% ct_exclude)), 0.9, 0.4))

df_metrics <- dplyr::select(df, cluster, cell_count, log10_cell_count, n_frag_per_cell, dataset, cluster_alpha) %>% distinct()

# plot
links_by_frag <- plot_scatter(df, "log10_unique_fragments", "n_enh_gene_links", label_key, dataset_key)
ggsave(file.path(out_dir, "multiome_and_atac_log_frag_by_links.pdf"), links_by_frag, width=4, height=4)

genes_by_frag <- plot_scatter(df, "log10_unique_fragments", "n_genes_with_enh", label_key, dataset_key)
ggsave(file.path(out_dir, "multiome_and_atac_frag_by_num_genes.pdf"), genes_by_frag, width=4, height=4)

size_by_num_enh <- plot_scatter(df, "n_enh_elements", "mean_enh_width", label_key, dataset_key)
ggsave(file.path(out_dir, "multiome_and_atac_enh_n_by_width.pdf"), size_by_num_enh, width=4, height=4)

num_enh_by_num_genes <- plot_scatter(df, "mean_enh_per_gene", "mean_genes_per_enh", label_key, dataset_key)
ggsave(file.path(out_dir, "multiome_and_atac_enh_per_gene_per_enh.pdf"), num_enh_by_num_genes, width=4, height=4)

cells_by_frag <- plot_scatter_simple(df_metrics, "log10_cell_count", "n_frag_per_cell", label_key)
ggsave(file.path(out_dir, "all_frag_by_log_cell.pdf"), cells_by_frag, width=4, height=4)

dist_to_tss <- plot_distribution(df, "mean_log10_dist_to_tss", label_key) 
ggsave(file.path(out_dir, "multiome_and_atac_dist_to_tss.pdf"), dist_to_tss, width=4, height=2.5)

gr =  plot_grid(cells_by_frag, links_by_frag, size_by_num_enh,
	dist_to_tss, num_enh_by_num_genes, genes_by_frag,
	nrow = 2, ncol = 3, rel_widths = c(1,1,1), rel_heights = c(1,1), align = "hv", axis = "lrtb")
ggsave2(file.path(out_dir, "multiome_and_atac_GRID.pdf"), gr, width=8, height=5)


# save stats file for supp table
if (TRUE) {
	stats <- fread(stats_file) %>% 
		pivot_longer(cols = -c(cluster, model_name), names_to = "metric", values_to = "value") %>% 
		mutate(model_name = ifelse(model_name == "multiome_powerlaw_v3", "scE2G_Multiome", "scE2G_ATAC"))
	
	fwrite(stats, file.path(out_dir, "properties_long_format.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
}