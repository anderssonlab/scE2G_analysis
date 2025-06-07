# libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(data.table)
  library(cowplot)
  library(ggdist)
})

add_stats_columns <- function(df, stats){
	# ratio	sample_id	num_umi	num_frag	frag_per_cell	umi_per_cell	num_cells
	df = left_join(df, stats, by="cluster") %>%
		rename(total_frag = num_frag, cell_count = num_cells, total_umi = num_umi,
			atac_frag_per_cell = frag_per_cell, rna_umi_per_cell = umi_per_cell)
	df$total_frag = df$cell_count * df$atac_frag_per_cell
	df$log10_total_frag = log10(df$total_frag)
	df$total_umi = df$cell_count * df$rna_umi_per_cell
	df$log10_total_umi = log10(df$total_umi)
	df$log10_cell_count = log10(df$cell_count)
	df$true_ratio = df$total_frag / df$total_umi

	return(df)
}


make_comparison_df <- function(df, stats, model1, model2){
    # read and process performance df... possible colnames:
		# df = pd.DataFrame(columns = ['cluster', 'model', 'AUPRC', 'AUPRC_95CI_low', 'AUPRC_95CI_high',
		# 						 'precision_70_pct_recall', 'precision_70_pct_recall_95CI_low', 'precision_70_pct_recall_95CI_high', 'threshold_70_pct_recall', 
		# 						 'precision_50_pct_recall', 'precision_50_pct_recall_95CI_low', 'precision_50_pct_recall_95CI_high', 'threshold_50_pct_recall', 
		# 						 'precision_model_threshold', 'precision_model_threshold_95CI_low', 'precision_model_threshold_95CI_high',
		# 						 'recall_model_threshold', 'recall_model_threshold_95CI_low', 'recall_model_threshold_95CI_high',
		# 						 'pct_missing_elements'])
	df = dplyr::select(df, model, pct_missing_total, AUPRC, precision_model_threshold, recall_model_threshold,
			total_frag, atac_frag_per_cell, cell_count) %>%
		dplyr::filter(model %in% c(model1, model2)) %>%
		mutate(model_id = ifelse(model==model1, "model1", "model2")) %>%
		dplyr::select(-c(model)) %>%
		pivot_wider(names_from=model_id, values_from=c(AUPRC, precision_model_threshold, recall_model_threshold, pct_missing_total)) %>%
		mutate(across(ends_with("_model1"), .names="{col}_diff") - across(ends_with("_model2"))) %>%
		rename_with(~ sub("_model1_diff", "", .), ends_with("_diff")) #%>%
		#add_stats_columns(stats)

		print(df)

    return(df)
}

plot_all_performance_metrics <-function(df_all, xvar, xlabel, cp, models_include, model_id, outdir){
	# process inputs
	df = dplyr::filter(df_all, model %in% models_include)
	df_dist = dplyr::filter(df_all, model == "distanceToTSS") %>% 
		select(model, model_name, AUPRC, precision_model_threshold, recall_model_threshold)
	df_abc = data.frame(list(model = "scATAC_ABC", model_name = "ABC (A=scATAC, C=power law)", AUPRC = 0.482, precision_model_threshold = 0.419, recall_model_threshold = 0.7))
	#df_abc = dplyr::filter(df_all, model == "scATAC_ABC", cluster == "atac_33387_rna_17961_cells_7821")
	df_lines = rbind(df_dist, df_abc)

	df <- df[sample(nrow(df)), ] # shuffle rows for when plotting >1 model

	# plotting params
	if (xvar %in% c("log10_total_frag")){
		xmin = min(df[[xvar]] - 0.05*min(df[[xvar]]))
		xmax = max(df[[xvar]] + 0.05*max(df[[xvar]]))
		frag_lines = c(2e6) %>% log10()
	} else if (xvar == "log10_total_umi") {
		xmin = min(df[[xvar]] - 0.05*min(df[[xvar]]))
		xmax = max(df[[xvar]] + 0.05*max(df[[xvar]]))
		frag_lines = c(1e6) %>% log10()
	} else if (xvar=="log10_cell_count") {
		xmin = min(df[[xvar]] - 0.05*min(df[[xvar]]))
		xmax = max(df[[xvar]])
		#frag_lines = c(100, 400, 1600, 6400) %>% log10()
		frag_lines=c(100) %>% log10()
	} else {
		xmin = 0
		xmax = max(df[[xvar]])
		frag_lines = c(0)
	}

	yvar = "AUPRC"
	ylabel = "AUPRC"
	auprc = ggplot(df, aes(x=df[[xvar]], y=!!sym(yvar), color=model_name)) +
		geom_vline(xintercept=frag_lines, linetype='dashed', color='#c5cad7') + # fragment count ref
		geom_hline(data=df_lines, aes(yintercept=!!sym(yvar), color=model_name), linetype="11") + # distance
		geom_point(alpha=0.7, shape = 16, size = 2) +
		labs(x = xlabel, y=ylabel) + ggtitle(ylabel) +
		xlim(c(xmin, xmax)) + ylim(c(0, 0.7)) +
		scale_color_manual(values=cp, name="Predictor") +
		theme_classic() + theme(axis.text = element_text(size = 7, color = "#000000"),
			axis.title = element_text(size = 8), axis.ticks = element_line(color = "#000000"), legend.position='None', aspect.ratio = 1)

	yvar = "precision_model_threshold"
	ylabel = "Precision at model threshold"
	prec = ggplot(df, aes(x=df[[xvar]], y=!!sym(yvar), color=model_name)) +
		geom_vline(xintercept=frag_lines, linetype='dashed', color='#c5cad7') + # fragment count ref
		geom_hline(data=df_lines, aes(yintercept=!!sym(yvar), color=model_name), linetype="11") + # distance
		geom_point(alpha=0.7, shape = 16, size = 2) +
		labs(x = xlabel, y=ylabel) + ggtitle("Precision") +
		xlim(c(xmin, xmax)) + ylim(c(0, 0.7)) +
		scale_color_manual(values=cp) +
		theme_classic() + theme(axis.text = element_text(size = 7, color = "#000000"),
			axis.title = element_text(size = 8), axis.ticks = element_line(color = "#000000"), legend.position='None', aspect.ratio = 1)

	yvar = "recall_model_threshold"
	ylabel = "Recall at model threshold"
	recall = ggplot(df, aes(x=df[[xvar]], y=!!sym(yvar), color=model_name)) +
		geom_vline(xintercept=frag_lines, linetype='dashed', color='#c5cad7') + # fragment count ref
		geom_hline(data=df_lines, aes(yintercept=!!sym(yvar), color=model_name), linetype="11") + # distance
		geom_point(alpha=0.7, shape = 16, size = 2) +
		labs(x = xlabel, y=ylabel) + ggtitle("Recall") +
		xlim(c(xmin, xmax)) + ylim(c(0, 0.75)) +
		scale_color_manual(values=cp) +
		theme_classic() + theme(axis.text = element_text(size = 7, color = "#000000"),
			axis.title = element_text(size = 8), axis.ticks = element_line(color = "#000000"), legend.position='None', aspect.ratio = 1)

	 # assemble and save
	outfile = file.path(outdir, paste0(model_id, "_performance_by_", xvar, ".pdf"))
	for_legend = auprc + theme(legend.position='top' ) +guides(color = guide_legend(nrow = 1))
	legend_h = as_grob(cowplot::get_plot_component(for_legend, "guide-box-top"))

	prow1 = plot_grid(auprc, prec, recall, nrow=1, rel_widths=c(1,1,1))
	grid = plot_grid(prow1, legend_h, nrow=2, rel_heights=c(1, 0.1))
	h = 2.75; w = 7

	ggsave2(outfile, plot=grid, width=w, height=h)
}

plot_heatmap_over_array <- function(df, fill_var, title, colors, diverge=FALSE, x_var = "atac_frag_per_cell", y_var = "cell_count") {
	if (fill_var %in% c("pct_missing", "total_frag", "total_umi")){
		df$fill_values = log10(df[[fill_var]])
		colors = c("#f6eff7", "#e5e5e9", "#c5cad7", "#96a0b3", "#6e788d", "#435369")
	} else if (!diverge) {
		df$fill_values = exp(df[[fill_var]])
	} else {
		df$fill_values =df[[fill_var]]
		df[[fill_var]] = round(df$fill_values, 2)
	}

	label_key <- c(atac_frag_per_cell = "Fragments per cell",
		cell_count = "# cells",
		rna_umi_per_cell = "UMIs per cell")
	
	 if (diverge) { # make color scale balanced
        #colors =  c("#af8dc3", "#f7f7f7", "#7fbf7b") # purple - white - green 792374", "#00647
		colors =  c("#006479", "#f7f7f7", "#792374") # teal - white - purple
		max_abs = max(0.04, df$fill_values)
        lims = c(-max_abs, max_abs)
    } else {
		lims = NULL
	}

    hm = ggplot(df, aes(x=!!sym(x_var), y=!!sym(y_var), fill=fill_values)) + 
        geom_tile() + 
        geom_text(aes(label =!!sym(fill_var)), size=2.5) +
        xlab(label_key[x_var]) + ylab(label_key[y_var]) + ggtitle(title) +
        coord_fixed(1) + 
        scale_fill_gradientn(colors=colors, oob=scales::squish, na.value="#FFFFFF", limits=lims) + 
        theme_classic() + theme(axis.text = element_text(size = 7, color = "#000000"), axis.title = element_text(size = 8), 
			axis.text.x = element_text(angle = 60, hjust=1), plot.title = element_text(size=10),
			axis.ticks = element_line(color = "#000000"), legend.position='None')

    return(hm)
}

plot_performance_heatmaps <- function(df, colors, output_dir, plot_id, comparison=FALSE, x_var = "atac_frag_per_cell", y_var = "cell_count") {
    # format fill vars
	df$log10_total_frag = round(log10(df$total_frag), 2)
    df$total_frag = round(df$total_frag/1e6, 1)
    df$pct_missing = round(df$pct_missing_total*100, 0)

	if (!comparison){
		df$AUPRC = round(df$AUPRC, 2)
		df$precision_model_threshold = round(df$precision_model_threshold, 2)
		df$recall_model_threshold= round(df$recall_model_threshold, 2)
	}

	df$atac_frag_per_cell = factor(df$atac_frag_per_cell, levels=unique(sort(df$atac_frag_per_cell)))
	df$cell_count = factor(df$cell_count, levels=unique(sort(df$cell_count)))

	if (x_var == "rna_umi_per_cell") {
		df$total_umi = round(df$total_umi/1e6, 1)
		df$rna_umi_per_cell = factor(df$rna_umi_per_cell, levels=unique(sort(df$rna_umi_per_cell)))

		total_frag = plot_heatmap_over_array(df, "total_umi", "Total RNA UMIs (M)", colors, FALSE, "rna_umi_per_cell")
		pct_missing = plot_heatmap_over_array(df, "pct_missing", "% CRISPR elements with no predictions", rev(colors), FALSE, "rna_umi_per_cell")
		auprc = plot_heatmap_over_array(df, "AUPRC", "AUPRC", colors, diverge=comparison, "rna_umi_per_cell")
		prec = plot_heatmap_over_array(df, "precision_model_threshold", "Precision at model threshold", colors, diverge=comparison, "rna_umi_per_cell")
		recall  = plot_heatmap_over_array(df, "recall_model_threshold", "Recall at model threshold", colors, diverge=comparison, "rna_umi_per_cell")
	}

	if (x_var == "atac_frag_per_cell") {
		total_frag = plot_heatmap_over_array(df, "total_frag", "Total ATAC fragments (M)", colors, FALSE)
		pct_missing = plot_heatmap_over_array(df, "pct_missing", "% CRISPR elements with no predictions", rev(colors), FALSE)
		auprc = plot_heatmap_over_array(df, "AUPRC", "AUPRC", colors, diverge=comparison)
		prec = plot_heatmap_over_array(df, "precision_model_threshold", "Precision at model threshold", colors, diverge=comparison)
		recall  = plot_heatmap_over_array(df, "recall_model_threshold", "Recall at model threshold", colors, diverge=comparison)
	}


    r1 = plot_grid(NULL, total_frag, pct_missing, NULL, nrow=1, rel_widths=c(1, 2, 2, 1)) 
    r2 = plot_grid(auprc, prec, recall, nrow=1) 

	if (comparison){
		output_file = file.path(output_dir, paste0(plot_id, "_", x_var, "_", y_var, "_comparison_heatmaps.pdf"))
		ggsave2(output_file, r2, width=10, height=3.25)
	} else {
		output_file = file.path(output_dir, paste0(plot_id,  "_", x_var, "_", y_var, "_performance_heatmaps.pdf"))
		grid = plot_grid(r1, r2, nrow=2, rel_heights=c(1,1), align = "hv", axis = "lrtb")
		ggsave2(output_file, grid, width=10, height=6.5)
	}
}

plot_downsample_dimensions <- function(df, x_var, y_var, out_dir) {
	label_key <- c(atac_frag_per_cell = "Fragments per cell",
		cell_count = "# cells",
		rna_umi_per_cell = "UMIs per cell",
		log10_total_frag = "log10 (total ATAC fragments)",
		log10_total_umi = "log10 (total RNA UMIs)")
	threshold_key <- c(atac_frag_per_cell = 0,
		cell_count = 100,
		rna_umi_per_cell = 0, 
		log10_total_frag = log10(2e6),
		log10_total_umi = log10(1e6))

	p <- ggplot(df, aes(x = !!sym(x_var), y = !!sym(y_var))) +
		geom_vline(xintercept=threshold_key[x_var], linetype='dashed', color='#c5cad7') + # fragment count ref
		geom_hline(yintercept=threshold_key[y_var], linetype='dashed', color='#c5cad7') + # fragment count ref
		geom_point(alpha=0.7, shape = 16, size = 2, color = "#435369") +
		labs(x = label_key[x_var], y = label_key[y_var]) +
		theme_classic() + theme(axis.text = element_text(size = 7, color = "#000000"),
			axis.title = element_text(size = 8), axis.ticks = element_line(color = "#000000"),
			legend.position='None', aspect.ratio = 1)

	out_file <- file.path(out_dir, paste0("downsample_dimsensions_", x_var, "_", y_var, ".pdf"))
	ggsave(out_file, p, width=3, height=3)
}

#################### RUN


## input files
# results_dirs <- c("2025_0209_downsample_ratio1", "2025_0209_downsample_ratio2", "2025_0209_downsample_ratio5", "2025_0209_downsample_ratio10")
# input_files <- file.path("/oak/stanford/groups/engreitz/Users/sheth/sc-E2G/results", results_dirs, "crispr_benchmarking_performance_summary.tsv")

## corrected multiome model downsample
results_dirs <- c("2025_0319_downsample_multiome_ratio1", "2025_0319_downsample_multiome_ratio2", "2025_0319_downsample_multiome_ratio5", "2025_0319_downsample_multiome_ratio10")
input_files <- file.path("/oak/stanford/groups/engreitz/Users/sheth/scE2G_temp/scE2G/results", results_dirs, "crispr_benchmarking_performance_summary.tsv")

# relevant columns: cluster model   AUPRC   AUPRC_95CI_low  AUPRC_95CI_high  precision_model_threshold       precision_model_threshold_95CI_low      precision_model_threshold_95CI_high     recall_model_threshold  recall_model_threshold_95CI_low recall_model_threshold_95CI_high        pct_missing_elements
stats_file = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0124_new_downsample/config/downsample_stats.tsv" # ratio	sample_id	num_umi	num_frag	frag_per_cell	umi_per_cell	num_cells
output_dir = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0124_new_downsample/apply_across_downsamples/"

dir.create(output_dir)

cp_df = data.frame(model = c("multiome_powerlaw_v3","scATAC_powerlaw_v3", "distanceToTSS", "multiome_powerlaw_v3_noTPMfilter", "scATAC_ABC"),
								hex = c("#792374", "#006479", "#6e788d",
									"#b778b3", "#0096a0"),
								model_name = c("scE2G (Multiome)", "scE2G (ATAC)", "Distance to TSS",
									"scE2G (Multome, no TPM filter)", "ABC (A=scATAC, C=power law)"))
cp = cp_df$hex; names(cp) = cp_df$model_name

## process input
stats = fread(stats_file, sep="\t") %>%
	rename(cluster = sample_id)

df_all = lapply(input_files, fread) %>%
	rbindlist() %>%
	as_tibble() %>%
	add_stats_columns(stats) %>% 
	left_join(cp_df, by="model")

## plot performance by fragment scatterplots
# plot by fragments, UMIs, cells, and ratio?
if (TRUE) {
	this_outdir <- file.path(output_dir, "performance_metrics_scatter")
	dir.create(this_outdir)
	x_vars <- c("log10_total_frag", "log10_total_umi", "log10_cell_count", "true_ratio",
		"atac_frag_per_cell", "rna_umi_per_cell")
	x_labels <- c("log10 (# unique ATAC fragments)", "log10 (# RNA UMIs)", "log10 (# cells)", "ATAC fragment : RNA UMI ratio",
		"Fragments per cell", "UMIs per cell")


	model_set1 <- c("multiome_powerlaw_v3", "scATAC_powerlaw_v3")
	#model_set2 <- c("multiome_powerlaw_v2", "multiome_powerlaw_v2_noTPMfilter")
	model_vars <- list(multiome_only = c("multiome_powerlaw_v3"))
	#model_vars <- list(multiome_versus_atac = model_set1, multiome_only=c("multiome_powerlaw_v3"), scATAC_only=c("scATAC_powerlaw_v3"))
		#multiome_TPM = model_set2, multiome_noTPMfilter_only=c("multiome_powerlaw_v3_noTPMfilter"), )
	#models <- c("multiome_powerlaw_v2", "scATAC_powerlaw_v2", "multiome_powerlaw_v2_noTPMfilter", "scATAC_ABC")

	for (m in seq_along(model_vars)) {
		models_include <- model_vars[[m]]
		model_id <- names(model_vars)[m]
		message(model_id)

		for (n in seq_along(x_vars)) {
			xvar <- x_vars[n]
			xlabel <- x_labels[n]
			message(xvar)

			plot_all_performance_metrics(df_all, xvar, xlabel, cp, models_include, model_id, this_outdir)
		}
	}
}


## heatmaps of performance per model
colors = c("#f6eff7","#bdc9e1", "#67a9cf","#1c9099", "#016c59") # white - teal (old performance colors)
purples <- c("#f6eff7", "#e9d3ea", "#d3a9ce", "#b778b3", "#a64791", "#792374") #, "#430b4e")
teals <- c("#f6eff7", "#cae5ee", "#96ced3", "#49bcbc", "#0096a0", "#006479") #, "#003648")

df_m = dplyr::filter(df_all, model=="multiome_powerlaw_v3")
#df_s = dplyr::filter(df_all, model=="scATAC_powerlaw_v3")

# df_s %>% arrange(AUPRC) %>% filter(total_frag > 2e6) %>% select(cluster, total_frag, total_umi, ratio) %>% head()
# df_m2 = dplyr::filter(df_all, model=="multiome_powerlaw_v2_noTPMfilter")

if (TRUE) {
	this_outdir <- file.path(output_dir, "performance_heatmaps")
	dir.create(this_outdir)
	ratios <- c(1, 2, 5, 10, "all")
	
	for (ratio_selected in ratios) {
		if (ratio_selected == "all") {
			df_m_this <- group_by(df_m, atac_frag_per_cell, cell_count) %>%
				summarize(AUPRC = mean(AUPRC), precision_model_threshold = mean(precision_model_threshold), recall_model_threshold = mean(recall_model_threshold),
					pct_missing_total = mean(pct_missing_total), total_frag = mean(total_frag))
			# df_s_this <- filter(df_s, cluster != "atac_33387_rna_16694_cells_500") %>% 
			# 	group_by(atac_frag_per_cell, cell_count) %>%
			# 	summarize(AUPRC = mean(AUPRC), precision_model_threshold = mean(precision_model_threshold), recall_model_threshold = mean(recall_model_threshold),
			# 		pct_missing_total = mean(pct_missing_total), total_frag = mean(total_frag))
		} else {
			df_m_this <- filter(df_m, ratio == ratio_selected)
			# df_s_this <- filter(df_s, ratio == ratio_selected)
		}

		plot_performance_heatmaps(df_m_this, purples, this_outdir, paste0("multiome_powerlaw_v3_", ratio_selected), comparison=FALSE)
		# plot_performance_heatmaps(df_s_this, teals, this_outdir, paste0("scATAC_powerlaw_v3_", ratio_selected), comparison=FALSE)
		#plot_performance_heatmaps(df_m2, colors, output_dir, "multiome_powerlaw_v2_noTPMfilter", comparison=FALSE)
	}
}

# heatmaps of UMI x cell perfomrance for multiome model
if (TRUE) {
	this_outdir <- file.path(output_dir, "performance_heatmaps")
	dir.create(this_outdir)
	ratios <- c(1, 2, 5, 10)
	
	for (ratio_selected in ratios) {
		df_m_this <- filter(df_m, ratio == ratio_selected)
		plot_performance_heatmaps(df_m_this, purples, this_outdir, paste0("multiome_powerlaw_v3_", ratio_selected),
			comparison=FALSE, x_var = "rna_umi_per_cell")
	}
}
	


## heatmaps of multiome versus scATAC performance

#df_comp2  = make_comparison_df(df_ratio, stats, model1="multiome_powerlaw_v2_noTPMfilter", model2="scATAC_powerlaw_v2")
#df_comp3  = make_comparison_df(df_ratio, stats, model1="multiome_powerlaw_v2", model2="multiome_powerlaw_v2_noTPMfilter")

if (FALSE) {
	this_outdir <- file.path(output_dir, "comparison_heatmaps")
	dir.create(this_outdir)
	ratios <- c(1, 2, 5, 10, "all")
	
	for (ratio_selected in ratios) {
		if (ratio_selected == "all") {
			df_ratio <-  filter(df_all, !(cluster == "atac_33387_rna_16694_cells_500" & model == "scATAC_powerlaw_v3")) %>%
			group_by(model, model_name, atac_frag_per_cell, cell_count) %>%
				summarize(AUPRC = mean(AUPRC), precision_model_threshold = mean(precision_model_threshold), recall_model_threshold = mean(recall_model_threshold),
					pct_missing_total = mean(pct_missing_total), total_frag = mean(total_frag)) %>%
			ungroup()
		} else {
			df_ratio <- dplyr::filter(df_all, ratio == ratio_selected)
		}

		df_comp = make_comparison_df(df_ratio, stats, model1="multiome_powerlaw_v3", model2="scATAC_powerlaw_v3")
		plot_performance_heatmaps(df_comp, colors, this_outdir, paste0("multiome_minus_scATAC_ratio_", ratio_selected), comparison=TRUE)

	}

	#plot_performance_heatmaps(df_comp2, colors, this_outdir, paste0("multiome_noTPM_minus_scATAC_ratio_", ratio_selected), comparison=TRUE)
	#plot_performance_heatmaps(df_comp3, colors, this_outdir, paste0("multiome_minus_multiome_noTPM_ratio_", ratio_selected), comparison=TRUE)
}

if (FALSE) {
	vars <- c("log10_total_frag", "log10_total_umi", "cell_count", "atac_frag_per_cell", "rna_umi_per_cell")
	df_dim <- df_all %>% select(all_of(vars)) %>% distinct()

	plot_downsample_dimensions(df_dim, "log10_total_frag", "log10_total_umi", output_dir)
	plot_downsample_dimensions(df_dim, "atac_frag_per_cell", "rna_umi_per_cell", output_dir)

}

