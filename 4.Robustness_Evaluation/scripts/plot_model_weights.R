# libraries
# libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(data.table)
  library(cowplot)
})

add_stats_columns <- function(df, stats){
	df = left_join(df, stats, by="cluster") 
	df$total_frag = df$cell_count * df$atac_frag_per_cell
	df$log10_total_frag = log10(df$total_frag)
	df$total_umi = df$cell_count * df$rna_umi_per_cell
	df$log10_total_umi = log10(df$total_umi)
	df$log10_cell_count = log10(df$cell_count)

	return(df)
}

get_model_weights <- function(stats_df, pred_dir, model_name, output_dir){
    df = stats_df
    df$model_weight_file = file.path(pred_dir, df$cluster, paste0(model_name, "_", df$cluster), "model", "model_coefficients.tsv")
    df_merge = dplyr::select(df, cluster, cell_count, atac_frag_per_cell) %>%
		mutate(total_frag = cell_count * atac_frag_per_cell)

    for (i in 1:nrow(df)){
        df_this = fread(df$model_weight_file[i])
        df_this = dplyr::filter(df_this, test_chr=="none") # select "full model" weights
        df_this$cluster = df$cluster[i]
        if (i==1) {df_mw = df_this} else (df_mw = rbind(df_mw, df_this))
    }

    df_mw = dplyr::left_join(df_mw, df_merge, by="cluster")
    write.table(df_mw, file.path(output_dir, paste0(model_name, "_", "model_weights.tsv")),  quote=FALSE, row.names=FALSE, sep="\t")
    return(df_mw)
}

plot_model_weights <- function(df_mw, col_var, legend_title, model_name, output_dir){
	#cols  = c("#f6eff7","#bdc9e1", "#67a9cf","#1c9099", "#016c59") # white - teal
	cols = c("#6e788d", "#435369", "#1c2a43") # greys
	atac_thresh = 2

    g = ggplot(df_mw, aes(x=total_frag/1E6, y=coefficient, color=df_mw[[col_var]])) +
        geom_point(size=1, alpha=0.8, shape = 16) +
		geom_vline(xintercept=atac_thresh, linetype="dashed", color="#96a0b3") +
        scale_x_log10() +
        facet_wrap(vars(feature), ncol=4, axes="all") + 
        xlab('Total ATAC fragments (M)') + ylab("Model weight") +
		scale_color_gradientn(colors=cols, na.value="#FFFFFF", name=legend_title) + 
        theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position='top', strip.background = element_blank(), aspect.ratio=1)

    ggsave(file.path(output_dir, paste0(model_name, "_", "weights_log_fill_", col_var, ".pdf")), g, width = 6, height = 4)
}

stats_file = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0524_stress_testing/downsample_stats.tsv"
pred_dir = "/oak/stanford/groups/engreitz/Users/sheth/sc-E2G/results/2024_0829_downsample"
output_dir = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0710_downsample_for_paper/train_across_downsamples"

dir.create(output_dir)

stats = fread(stats_file, sep="\t") %>%
	dplyr::filter(!(cluster %in% c("K562_IGVF", "K562_Wang")))

arc = get_model_weights(stats, pred_dir, "arc", output_dir)
abck = get_model_weights(stats, pred_dir, "abck", output_dir)
atac = get_model_weights(stats, pred_dir, "atac", output_dir)

plot_model_weights(arc, "cell_count", "# cells", "arc", output_dir)
plot_model_weights(abck, "cell_count", "# cells", "abck", output_dir)
plot_model_weights(atac, "cell_count", "# cells", "atac", output_dir)

#plot_model_weights(df_mw, "frag_per_cell", "Greens", "Fragments per cell", output_dir)
