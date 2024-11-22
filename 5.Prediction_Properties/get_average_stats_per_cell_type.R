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
	model_cp <- c(multiome_powerlaw_v2 = "#792374", scATAC_powerlaw_v2 = "#006479")

	# scatter plot
	s <- ggplot(df, aes(x = !!sym(x_metric), y = !!sym(y_metric), alpha = cluster_alpha)) +
		geom_vline(xintercept = vline_loc, alpha = vline_alpha, linetype = "dashed", color = "#96a0b3") +
		geom_point(aes(color = model_name), size = 1.25) +
		labs(x = x_label, y = label_key[y_metric]) +
		scale_color_manual(values = model_cp) +
		scale_alpha_identity() +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), aspect.ratio = 1,
				legend.position = "None")

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
	model_cp <- c(multiome_powerlaw_v2 = "#792374", scATAC_powerlaw_v2 = "#006479")
	x_label <- paste0(label_key[metric], "\n", n_label)
	
	x <- ggplot(df, aes(x = !!sym(metric), y = model_name, fill = model_name, alpha = cluster_alpha)) +
		geom_dots(side = "top", layout = "swarm", linewidth = 0) + 
		labs(x = x_label, y = NULL) +
		scale_fill_manual(values = model_cp) +
		scale_alpha_identity() +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), axis.text.y = element_blank(), legend.position = "None")
}

plot_scatter_simple <- function(df, x_metric, y_metric, label_key) {
	# format data
	n_label <- paste0("N = ", nrow(df))
	x_label <- paste0(label_key[x_metric], "\n", n_label)



	# scatter plot
	s <- ggplot(df, aes(x = !!sym(x_metric), y = !!sym(y_metric), alpha = cluster_alpha)) +
		geom_point(aes(color = model_name), size = 1.25, color = "#1c2a43") +
		labs(x = x_label, y = label_key[y_metric]) + scale_alpha_identity() +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), aspect.ratio = 1,
				legend.position = "None")
}

##### MAIN

## get inputs from snakemake
# sample_key <- fread(snakemake@input$sample_key)
# score_threshold <- snakemake@params$score_threshold %>% as.numeric()
# score_column <- snakemake@params$score_column
# include_promoters <- ifelse(snakemake@params$include_promoters %in% c("TRUE", "True"), TRUE, FALSE) 
stats_file <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0916_global_properties/results_v1/multiome_pred_stats.tsv"
out_dir <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0916_global_properties/results_v1"
model_name <- "multiome_powerlaw_v2"

# labels 
metric_names <- c("num_unique_fragments_M", "num_cells", "num_enh", "num_genes_with_enh", "num_enh_gene_links",
		"mean_num_genes_per_enh", "mean_num_enh_per_gene", "mean_log10_dist_to_tss", "mean_enh_region_size",
		"num_genes_considered", "num_genes_not_expressed", "log10_unique_fragments", "num_frag_per_cell", "log10_num_cells")
metric_labels <- c("# unique ATAC fragments (M)", "# cells", "# enhancer elements", "# genes with 1+ enhancer", "# enhancer-gene links",
	"Mean # genes per enhancer", "Mean # enhancers per gene", "Mean log10 distance to TSS (bp)", "Mean enhancer width (bp)",
	"# genes considered", "# genes not expressed", "log10(# unique ATAC fragments)", "Mean # ATAC fragments per cell", "log10(# cells)")
label_key <- setNames(metric_labels, metric_names) # name: label

# biosample color key
datasets <- c("K562", "GM12878", "BMMC_all", "PBMC_all", "Islets_all", "BMMC5", "PBMC5", "Islets", "BMMC22", "PBMC9")
dataset_cols <- c(rep("#1c2a43", 5), rep("#435369", 3), rep("#6e788d", 2))
dataset_key <- setNames(dataset_cols, datasets)

# read in and format data
sample_key <- fread("/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0916_global_properties/config/pred_sample_key.tsv") %>%
	dplyr::select(biosample, biosample_name, dataset)

stats <- fread(stats_file) %>% 
	dplyr::filter(model_name == model_name) %>%
	left_join(sample_key, by=c("cell_cluster" = "biosample"))

# add atac data
atac_file <-  "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0916_global_properties/results_v1/atac_pred_stats.tsv"
atac <- fread(atac_file) %>% dplyr::filter(model_name == "scATAC_powerlaw_v2") %>% left_join(sample_key, by=c("cell_cluster" = "biosample"))
stats <- rbind(stats, atac)

# filter to dataset
ds_include <- c("BMMC22", "PBMC9", "Islets")
stats_for_summ <- stats %>%
	dplyr::filter(dataset %in% ds_include, cell_cluster  != "BMMC22_ID2_hi_myeloid_prog", model_name == "multiome_powerlaw_v2") 

calcs <- pivot_wider(stats_for_summ, names_from = metric, values_from = value) %>%
	mutate(total_enh_size = num_enh * mean_enh_region_size,
		total_enh_log10_dist = num_enh_gene_links * mean_log10_dist_to_tss)
mean_size <- sum(calcs$total_enh_size/sum(calcs$num_enh))
mean_dist <- sum(calcs$total_enh_log10_dist/sum(calcs$num_enh_gene_links))
message("mean enhancer size: ", mean_size)
message("total number of enhancers: ", sum(calcs$num_enh))
message("mean enhancer log10 distance to tss: ", mean_dist)
message("total number of links: ", sum(calcs$num_enh_gene_links))

summary <- stats_for_summ %>%
	group_by(metric) %>%
	summarize(mean = mean(value), median = median(value), min = min(value), max = max(value), n_clusters = n())
fwrite(summary, file.path(out_dir, "multiome_stats_summary_granular_2M_noKG.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# pivot wider
df <- pivot_wider(stats, names_from = metric, values_from = value) %>%
	mutate(num_unique_fragments_M = num_unique_fragments/1e6,
		log10_unique_fragments = log10(num_unique_fragments),
		log10_num_cells = log10(num_cells),
		num_frag_per_cell = num_unique_fragments/num_cells,
		cluster_alpha= ifelse((dataset %in% ds_include & cell_cluster != "BMMC22_ID2_hi_myeloid_prog"), 0.9, 0.4))

df_metrics <- dplyr::select(df, cell_cluster, num_cells, log10_num_cells, num_frag_per_cell, dataset, cluster_alpha) %>% distinct()

# plot
links_by_frag <- plot_scatter(df, "log10_unique_fragments", "num_enh_gene_links", label_key, dataset_key)
ggsave(file.path(out_dir, "multiome_and_atac_log_frag_by_links.pdf"), links_by_frag, width=4, height=4)

genes_by_frag <- plot_scatter(df, "log10_unique_fragments", "num_genes_with_enh", label_key, dataset_key)
ggsave(file.path(out_dir, "multiome_and_atac_frag_by_num_genes.pdf"), genes_by_frag, width=4, height=4)

size_by_num_enh <- plot_scatter(df, "num_enh", "mean_enh_region_size", label_key, dataset_key)
ggsave(file.path(out_dir, "multiome_and_atac_enh_n_by_width.pdf"), size_by_num_enh, width=4, height=4)

num_enh_by_num_genes <- plot_scatter(df, "mean_num_enh_per_gene", "mean_num_genes_per_enh", label_key, dataset_key)
ggsave(file.path(out_dir, "multiome_and_atac_enh_per_gene_per_enh.pdf"), num_enh_by_num_genes, width=4, height=4)

cells_by_frag <- plot_scatter_simple(df_metrics, "log10_num_cells", "num_frag_per_cell", label_key)
ggsave(file.path(out_dir, "all_frag_by_log_cell.pdf"), cells_by_frag, width=4, height=4)

dist_to_tss <- plot_distribution(df, "mean_log10_dist_to_tss", label_key) 
ggsave(file.path(out_dir, "multiome_and_atac_dist_to_tss.pdf"), dist_to_tss, width=4, height=2.5)

gr =  plot_grid(cells_by_frag, links_by_frag, size_by_num_enh,
	dist_to_tss, num_enh_by_num_genes, genes_by_frag,
	nrow = 2, ncol = 3, rel_widths = c(1,1,1), rel_heights = c(1,1), align = "hv", axis = "lrtb")
ggsave2(file.path(out_dir, "multiome_and_atac_GRID.pdf"), gr, width=8, height=5)

