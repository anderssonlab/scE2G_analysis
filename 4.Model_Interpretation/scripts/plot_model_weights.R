# libraries
# libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(data.table)
  library(cowplot)
})

read_stats_df <- function(stats_fie){
	# ratio	sample_id	num_umi	num_frag	frag_per_cell	umi_per_cell	num_cells
	df <- fread(stats_file) %>% 
		rename(total_frag = num_frag, cell_count = num_cells, total_umi = num_umi,
			atac_frag_per_cell = frag_per_cell, rna_umi_per_cell = umi_per_cell,
            cluster = sample_id)
	df$total_frag = df$cell_count * df$atac_frag_per_cell
	df$log10_total_frag = log10(df$total_frag)
	df$total_umi = df$cell_count * df$rna_umi_per_cell
	df$log10_total_umi = log10(df$total_umi)
	df$log10_cell_count = log10(df$cell_count)
	df$true_ratio = df$total_frag / df$total_umi

	return(df)
}

read_one_file <- function(i, df) {
    df_this <- fread(df$model_weight_file[i]) %>% 
        dplyr::filter(test_chr == "none") %>% 
        mutate(cluster = df$cluster[i])
    return(df_this)
}

get_model_weights <- function(stats_df, results_dirs, model_name, output_dir){
    df <- stats_df %>% 
        mutate(results_dir = ifelse(ratio == 2, results_dirs[1], results_dirs[2]),
            model_id = paste0(model_name, "_", cluster),
            model_weight_file = file.path(results_dir, cluster, model_id, "model", "model_coefficients.tsv"))
    print(head(df))

    df_merge = dplyr::select(df, cluster, cell_count, atac_frag_per_cell, rna_umi_per_cell,
        log10_total_frag, log10_total_umi)

    df_mw <- lapply(1:nrow(df), read_one_file, df) %>% 
        rbindlist() %>% as_tibble() %>% 
        left_join(df_merge, by = "cluster")

    write.table(df_mw, file.path(output_dir, paste0(model_name, "_", "model_weights.tsv")),  quote=FALSE, row.names=FALSE, sep="\t")
    return(df_mw)
}

plot_model_weights <- function(df_mw, feature_table_file, model_name, output_dir){
	#cols  = c("#f6eff7","#bdc9e1", "#67a9cf","#1c9099", "#016c59") # white - teal
	greys <- c("f6eff7", "#e5e5e9", "#c5cad7", "#96a0b3", "#6e788d", "#435369", "#1c2a43")
    purples <- c("#e9d3ea", "#d3a9ce", "#b778b3", "#a64791", "#792374", "#430b4e")
    teals <- c("#cae5ee", "#96ced3", "#49bcbc", "#0096a0", "#006479", "#003648")
    greens <- c("#d7e5c5", "#a1ca78", "#5eb342", "#429130", "#1c6e2b", "#0e3716")
    blues <- c("#c5e5fb", "#9bcae9", "#5496ce", "#006eae", "#00488d", "#002359")

    colors <- case_when(model_name == "scATAC" ~ teals, model_name == "ARC" ~ purples, model_name == "ABCK" ~ blues)

	atac_thresh <- log10(2e6)
    rna_thresh <- log10(1e6)
    ft <- fread(feature_table_file)

    features_highlight <- c("ABC.Score", "Kendall", "ARC.E2G.Score")
    others <- setdiff(ft$feature, features_highlight) %>% sort()
    full_order <- c(features_highlight, others)

    df_mw <- df_mw %>% 
        left_join(ft, by = "feature") %>%
        mutate(color = ifelse(feature %in% features_highlight, colors[6], colors[3]),
            feature = factor(feature, levels = full_order, ordered = TRUE)) %>%
        arrange(feature)
    df_mw$nice_name <- factor(df_mw$nice_name, levels = unique(df_mw$nice_name), ordered = TRUE)

    atac = ggplot(df_mw, aes(x=log10_total_frag, y=coefficient)) +
        geom_point(aes(color = color), size=1, alpha=0.7, shape = 16) +
		geom_vline(xintercept=atac_thresh, linetype="dashed", color="#c5cad7") +
        facet_wrap(vars(nice_name), nrow=1, axes="all") + 
        xlab('log10 (# unique ATAC fragments)') + ylab("Model weight") +
        scale_color_identity() +
		#scale_color_gradientn(colors=cols, na.value="#FFFFFF", name=legend_title) + 
        theme_classic() + theme(axis.text = element_text(size = 7, color = "#000000"), axis.title = element_text(size = 8), #legend.position='top',
            axis.ticks = element_line(color = "#000000"), strip.background = element_blank(), aspect.ratio=1,
            legend.position = "none")

    rna = ggplot(df_mw, aes(x=log10_total_umi, y=coefficient)) +
        geom_point(aes(color = color), size=1, alpha=0.7, shape = 16) +
		geom_vline(xintercept=rna_thresh, linetype="dashed", color="#c5cad7") +
        facet_wrap(vars(nice_name), nrow=1, axes="all") + 
        xlab('log10 (# RNA UMIs)') + ylab("Model weight") +
        scale_color_identity() +
		#scale_color_gradientn(colors=cols, na.value="#FFFFFF", name=legend_title) + 
        theme_classic() + theme(axis.text = element_text(size = 7, color = "#000000"), axis.title = element_text(size = 8), #legend.position='top',
            axis.ticks = element_line(color = "#000000"), strip.background = element_blank(), aspect.ratio=1,
            legend.position = "none")
    
    g <- plot_grid(atac, rna, nrow = 2, align = "hv", rel_heights = c(1, 1))
    ggsave2(file.path(output_dir, paste0(model_name, "_", "model_weights.pdf")), g, width = 8, height = 4)

    return(list(atac, rna))
}

stats_file = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0124_new_downsample/config/downsample_stats.tsv" # ratio	sample_id	num_umi	num_frag	frag_per_cell	umi_per_cell	num_cells

results_names <- c("2025_0319_downsample_multiome_ratio2", "2025_0319_downsample_multiome_ratio10")
results_dirs <- file.path("/oak/stanford/groups/engreitz/Users/sheth/scE2G_temp/scE2G/results", results_names)
output_dir = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0124_new_downsample/train_across_downsamples"
dir.create(output_dir, showWarnings = FALSE)

abck_feature_file <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0626_feature_analysis_for_paper/feature_table_all_multiome_abck.tsv"
arc_feature_file <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_temp/scE2G/models/multiome_powerlaw_v3/feature_table.tsv"
atac_feature_file <- "/oak/stanford/groups/engreitz/Users/sheth/sc-E2G/models/scATAC_powerlaw_v3/feature_table.tsv"


stats = read_stats_df(stats_file) %>% 
    filter(ratio %in% c(2, 10))

# arc = get_model_weights(stats, results_dirs, "ARC", output_dir)
# abck = get_model_weights(stats, results_dirs, "ABCK", output_dir)
# atac = get_model_weights(stats, results_dirs, "scATAC", output_dir)

#models <- c("ARC", "ABCK", "scATAC")
models <- c("ARC", "ABCK")
ft_files <- c(arc_feature_file, abck_feature_file, atac_feature_file)
pl <- vector("list", length(models))
for (i in seq_along(models)) {
    model_name = models[i]
    mw_file <- file.path(output_dir, paste0(model_name, "_", "model_weights.tsv"))
    mw <- fread(mw_file)
    pl[[model_name]] <- plot_model_weights(mw, ft_files[i], model_name, output_dir)
}

plots <- pl[c("ABCK", "ARC")] %>% unlist(recursive = FALSE)

new_grid <- plot_grid(plotlist = plots, nrow = 4, align = "v", axis = "lrtb", rel_heights = c(1, 1, 1, 1))
ggsave2(file.path(output_dir, paste0("ABCK_ARC_", "_", "model_weights.pdf")), new_grid, width = 9, height = 6)


#plot_model_weights(df_mw, "frag_per_cell", "Greens", "Fragments per cell", output_dir)
