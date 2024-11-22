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
	df = left_join(df, stats, by="cluster") # model = multiome_7features or scATAC_6features
	df$total_frag = df$cell_count * df$atac_frag_per_cell
	df$log10_total_frag = log10(df$total_frag)
	df$total_umi = df$cell_count * df$rna_umi_per_cell
	df$log10_total_umi = log10(df$total_umi)
	df$log10_cell_count = log10(df$cell_count)

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
	df = dplyr::select(df, cluster, model, pct_missing_total, AUPRC, precision_model_threshold, recall_model_threshold) %>%
		dplyr::filter(model %in% c(model1, model2)) %>%
		mutate(model_id = ifelse(model==model1, "model1", "model2")) %>%
		dplyr::select(-c(model)) %>%
		pivot_wider(names_from=model_id, values_from=c(AUPRC, precision_model_threshold, recall_model_threshold, pct_missing_total)) %>%
		mutate(across(ends_with("_model1"), .names="{col}_diff") - across(ends_with("_model2"))) %>%
		rename_with(~ sub("_model1_diff", "", .), ends_with("_diff")) %>%
		add_stats_columns(stats)

		print(df)

    return(df)
}

plot_all_performance_metrics <-function(df_all, xvar, xlabel, cp, model_id, outdir){
	# process inputs
	df = dplyr::filter(df_all, model==model_id, dataset=="Xu")
	df_dist = dplyr::filter(df_all, model == "distanceToTSS")
	df_abc = dplyr::filter(df_all, model == "scATAC_ABC", cluster == "all_cells_all_atac")
	df_lines = rbind(df_dist, df_abc)

	# plotting params
	if (xvar %in% c("log10_total_frag", "log10_total_umi")){
		xmin = min(df[[xvar]] - 0.05*min(df[[xvar]]))
		xmax = max(df[[xvar]] + 0.05*max(df[[xvar]]))
		frag_lines = c(2e6) %>% log10()
	} else if (xvar=="log10_cell_count") {
		xmin = min(df[[xvar]] - 0.05*min(df[[xvar]]))
		xmax = max(df[[xvar]])
		#frag_lines = c(100, 400, 1600, 6400) %>% log10()
		frag_lines=c(0)
	}

	yvar = "AUPRC"
	ylabel = "AUPRC"
	auprc = ggplot(df, aes(x=df[[xvar]], y=!!sym(yvar), color=model_name)) +
		geom_vline(xintercept=frag_lines, linetype='dashed', color='#c5cad7') + # fragment count ref
		geom_hline(data=df_lines, aes(yintercept=!!sym(yvar), color=model_name), linetype="11") + # distance
		geom_point(alpha=0.8, shape = 16, size = 2) +
		labs(x = xlabel, y=ylabel) + ggtitle(ylabel) +
		xlim(c(xmin, xmax)) + ylim(c(0, 0.7)) +
		scale_color_manual(values=cp, name="Predictor") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position='None', aspect.ratio = 1)

	yvar = "precision_model_threshold"
	ylabel = "Precision at model threshold"
	prec = ggplot(df, aes(x=df[[xvar]], y=!!sym(yvar), color=model_name)) +
		geom_vline(xintercept=frag_lines, linetype='dashed', color='#c5cad7') + # fragment count ref
		geom_hline(data=df_lines, aes(yintercept=!!sym(yvar), color=model_name), linetype="11") + # distance
		geom_point(alpha=0.8, shape = 16, size = 2) +
		labs(x = xlabel, y=ylabel) + ggtitle("Precision") +
		xlim(c(xmin, xmax)) + ylim(c(0, 0.7)) +
		scale_color_manual(values=cp) +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position="None", aspect.ratio = 1)

	yvar = "recall_model_threshold"
	ylabel = "Recall at model threshold"
	recall = ggplot(df, aes(x=df[[xvar]], y=!!sym(yvar), color=model_name)) +
		geom_vline(xintercept=frag_lines, linetype='dashed', color='#c5cad7') + # fragment count ref
		#geom_hline(data=df_dist, aes(yintercept=!!sym(yvar), color=model_name), linetype="solid") + # distance
		geom_point(alpha=0.8, shape = 16, size = 2) +
		labs(x = xlabel, y=ylabel) + ggtitle("Recall") +
		xlim(c(xmin, xmax)) + ylim(c(0, 0.75)) +
		scale_color_manual(values=cp) +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position="None", aspect.ratio = 1)

	 # assemble and save
	outfile = file.path(outdir, paste0(model_id, "_performance_by_", xvar, "_v3.pdf"))
	for_legend = auprc + theme(legend.position='top' ) +guides(color = guide_legend(nrow = 1))
	legend_h = as_grob(cowplot::get_plot_component(for_legend, "guide-box-top"))

	prow1 = plot_grid(auprc, prec, recall, nrow=1, rel_widths=c(1,1,1))
	grid = plot_grid(prow1, legend_h, nrow=2, rel_heights=c(1, 0.1))
	h = 2.75; w = 7

	ggsave2(outfile, plot=grid, width=w, height=h)
}

plot_heatmap_over_array <- function(df, fill_var, title, colors, diverge=FALSE) {
	if (fill_var %in% c("pct_missing", "total_frag")){
		df$fill_values = log10(df[[fill_var]])
	} else if (!diverge) {
		df$fill_values = exp(df[[fill_var]])
	} else {
		df$fill_values =df[[fill_var]]
		df[[fill_var]] = round(df$fill_values, 2)
	}
	
	 if (diverge) { # make color scale balanced
        #colors =  c("#af8dc3", "#f7f7f7", "#7fbf7b") # purple - white - green 792374", "#00647
		colors =  c("#006479", "#f7f7f7", "#792374") # teal - white - magenta
		max_abs = max(0.04, df$fill_values)
        lims = c(-max_abs, max_abs)
    } else {
		lims = NULL
	}

    hm = ggplot(df, aes(x=atac_frag_per_cell, y=cell_count, fill=fill_values)) + 
        geom_tile() + 
        geom_text(aes(label =!!sym(fill_var)), size=2.5) +
        xlab('Fragments per cell') + ylab('# cells') + ggtitle(title) +
        coord_fixed(1) + 
        scale_fill_gradientn(colors=colors, oob=scales::squish, na.value="#FFFFFF", limits=lims) + 
        theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), 
													axis.text.x = element_text(angle = 60, hjust=1),
													plot.title = element_text(size=10), legend.position='None')

    return(hm)
}

plot_performance_heatmaps <- function(df, colors, output_dir, plot_id, comparison=FALSE) {
    # format fill vars
	df$log10_total_frag = round(log10(df$total_frag), 2)
    df$total_frag = round(df$total_frag/1e6, 1)
	df$total_umi = round(df$total_umi/1e6, 1)
    df$pct_missing = round(df$pct_missing_total*100, 0)

	if (!comparison){
		df$AUPRC = round(df$AUPRC, 2)
		df$precision_model_threshold = round(df$precision_model_threshold, 2)
		df$recall_model_threshold= round(df$recall_model_threshold, 2)
	}


	df$atac_frag_per_cell = factor(df$atac_frag_per_cell, levels=unique(sort(df$atac_frag_per_cell)))
	df$cell_count = factor(df$cell_count, levels=unique(sort(df$cell_count)))

    total_frag = plot_heatmap_over_array(df, "total_frag", "Total ATAC fragments (M)", colors, FALSE)
	#total_frag = plot_heatmap_over_array(df, "log10_total_frag", "log10 (total ATAC fragments)", colors, FALSE)
    pct_missing = plot_heatmap_over_array(df, "pct_missing", "% CRISPR elements with no predictions", rev(colors), FALSE)
    auprc = plot_heatmap_over_array(df, "AUPRC", "AUPRC", colors, diverge=comparison)
    prec = plot_heatmap_over_array(df, "precision_model_threshold", "Precision at model threshold", colors, diverge=comparison)
    recall  = plot_heatmap_over_array(df, "recall_model_threshold", "Recall at model threshold", colors, diverge=comparison)

    r1 = plot_grid(NULL, total_frag, pct_missing, NULL, nrow=1, rel_widths=c(1, 2, 2, 1)) 
    r2 = plot_grid(auprc, prec, recall, nrow=1) 

	if (comparison){
		output_file = file.path(output_dir, paste0(plot_id, "_comparison_heatmaps.pdf"))
		ggsave2(output_file, r2, width=10, height=3.25)
	} else {
		output_file = file.path(output_dir, paste0(plot_id, "_performance_heatmaps.pdf"))
		grid = plot_grid(r1, r2, nrow=2, rel_heights=c(1,1), align = "hv", axis = "lrtb")
		ggsave2(output_file, grid, width=10, height=6.5)
	}
}

plot_jaccard_similarity <- function(jacc, colors, output_dir) {
	## format data
	num_el = max(jacc$n_intersections_51)
	jacc$atac_frag_per_cell = factor(jacc$atac_frag_per_cell, levels=unique(sort(jacc$atac_frag_per_cell)))
	jacc$cell_count = factor(jacc$cell_count, levels=unique(sort(jacc$cell_count)))
	lims=c(0,1)

	jacc = jacc %>%
		mutate(jaccard_51_lab = round(jaccard_51, 2),
			jaccard_90_lab = round(jaccard_90, 2),
			frac_int_51 = n_intersections_51/num_el, frac_int_51_lab = round(frac_int_51, 2),
			frac_int_90 = n_intersections_90/num_el, frac_int_90_lab = round(frac_int_90, 2))

	j5 = ggplot(jacc, aes(x=atac_frag_per_cell, y=cell_count, fill=jaccard_51)) + 
        geom_tile() + 
        geom_text(aes(label =jaccard_51_lab), size=2.5) +
        xlab('Fragments per cell') + ylab('# cells') + ggtitle("Jaccard similarity (50% overlap)") +
        coord_fixed(1) + 
        scale_fill_gradientn(colors=colors, oob=scales::squish, na.value="#FFFFFF", limits=lims) + 
        theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), 
													axis.text.x = element_text(angle = 60, hjust=1),
													plot.title = element_text(size=10), legend.position='None')
		
	f5 = ggplot(jacc, aes(x=atac_frag_per_cell, y=cell_count, fill=frac_int_51)) + 
        geom_tile() + 
        geom_text(aes(label =frac_int_51_lab), size=2.5) +
        xlab('Fragments per cell') + ylab('# cells') + ggtitle("Fraction overlap (min. 50%)") +
        coord_fixed(1) + 
        scale_fill_gradientn(colors=colors, oob=scales::squish, na.value="#FFFFFF", limits=lims) + 
        theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), 
													axis.text.x = element_text(angle = 60, hjust=1),
													plot.title = element_text(size=10), legend.position='None')
													
	j9 = ggplot(jacc, aes(x=atac_frag_per_cell, y=cell_count, fill=jaccard_90)) + 
        geom_tile() + 
        geom_text(aes(label =jaccard_90_lab), size=2.5) +
        xlab('Fragments per cell') + ylab('# cells') + ggtitle("Jaccard similarity (90% overlap)") +
        coord_fixed(1) + 
        scale_fill_gradientn(colors=colors, oob=scales::squish, na.value="#FFFFFF", limits=lims) + 
        theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), 
													axis.text.x = element_text(angle = 60, hjust=1),
													plot.title = element_text(size=10), legend.position='None')

	f9 = ggplot(jacc, aes(x=atac_frag_per_cell, y=cell_count, fill=frac_int_90)) + 
        geom_tile() + 
        geom_text(aes(label =frac_int_90_lab), size=2.5) +
        xlab('Fragments per cell') + ylab('# cells') + ggtitle("Fraction overlap (min. 90%)") +
        coord_fixed(1) + 
        scale_fill_gradientn(colors=colors, oob=scales::squish, na.value="#FFFFFF", limits=lims) + 
        theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), 
													axis.text.x = element_text(angle = 60, hjust=1),
													plot.title = element_text(size=10), legend.position='None')

	output_file = file.path(output_dir, "candidate_element_overlap_heatmaps.pdf")
	row = plot_grid(j5, f5, j9, f9, nrow=1, rel_widths=c(1,1,1,1))
	ggsave2(output_file, row, width=10, height=3.25)
}


# input files
input_file  = "/oak/stanford/groups/engreitz/Users/sheth/sc-E2G/results/2024_0829_downsample/crispr_benchmarking_performance_summary.tsv"
# relevant columns: cluster model   AUPRC   AUPRC_95CI_low  AUPRC_95CI_high  precision_model_threshold       precision_model_threshold_95CI_low      precision_model_threshold_95CI_high     recall_model_threshold  recall_model_threshold_95CI_low recall_model_threshold_95CI_high        pct_missing_elements
stats_file = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0524_stress_testing/downsample_stats.tsv" # cluster	cell_count	rna_umi_per_cell	atac_frag_per_cell	dataset
output_dir = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0710_downsample_for_paper/v2_apply_across_downsamples/"
jaccard_file = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0710_downsample_for_paper/v2_apply_across_downsamples/jaccard.tsv"

dir.create(output_dir)

cp_df = data.frame(model = c("multiome_powerlaw_v2", "scATAC_powerlaw_v2", "distanceToTSS",
	"multiome_powerlaw_v2_noTPMfilter", "scATAC_ABC"),
								hex = c("#792374", "#006479", "#6e788d",
									"#b778b3", "#0096a0"),
								model_name = c("scE2G (Multiome)", "scE2G (ATAC)", "Distance to TSS",
									"scE2G (Multome, no TPM filter)", "ABC (A=scATAC, C=power law)"))
cp = cp_df$hex; names(cp) = cp_df$model_name

## process input
stats = fread(stats_file, sep="\t")
df_all = fread(input_file, sep="\t") %>%
	add_stats_columns(stats) %>% 
	left_join(cp_df, by="model")

jacc = fread(jaccard_file, sep="\t") %>% 
	mutate(cluster=pred_id) %>%
	add_stats_columns(stats)

## plot performance by fragment scatterplots
if (TRUE) {
	plot_all_performance_metrics(df_all, xvar="log10_total_frag", xlabel="log10 (total ATAC fragments)", cp = cp, model_id ="multiome_powerlaw_v2", outdir= output_dir)
	plot_all_performance_metrics(df_all, xvar="log10_total_frag", xlabel="log10 (total ATAC fragments)", cp = cp, model_id="scATAC_powerlaw_v2", outdir= output_dir)
	plot_all_performance_metrics(df_all, xvar="log10_total_frag", xlabel="log10 (total ATAC fragments)", cp = cp, model_id ="multiome_powerlaw_v2_noTPMfilter", outdir= output_dir)
	plot_all_performance_metrics(df_all, xvar="log10_cell_count", xlabel="log10 (# cells)", cp = cp, model_id ="multiome_powerlaw_v2", outdir= output_dir)
	plot_all_performance_metrics(df_all, xvar="log10_cell_count", xlabel="log10 (# cells)", cp = cp, model_id ="scATAC_powerlaw_v2", outdir= output_dir)
	plot_all_performance_metrics(df_all, xvar="log10_cell_count", xlabel="log10 (# cells)", cp = cp, model_id ="multiome_powerlaw_v2_noTPMfilter", outdir= output_dir)
	plot_all_performance_metrics(df_all, xvar="log10_total_umi", xlabel="log10 (total RNA UMIs)", cp = cp, model_id ="multiome_powerlaw_v2", outdir= output_dir)
	plot_all_performance_metrics(df_all, xvar="log10_total_umi", xlabel="log10 (total RNA UMIs)", cp = cp, model_id ="scATAC_powerlaw_v2", outdir= output_dir)
	plot_all_performance_metrics(df_all, xvar="log10_total_umi", xlabel="log10 (total RNA UMIs)", cp = cp, model_id ="multiome_powerlaw_v2_noTPMfilter", outdir= output_dir)
}

## heatmaps of performance per model
colors = c("#f6eff7","#bdc9e1", "#67a9cf","#1c9099", "#016c59") # white - teal
df_m = dplyr::filter(df_all, model=="multiome_powerlaw_v2", dataset=="Xu")
df_s = dplyr::filter(df_all, model=="scATAC_powerlaw_v2", dataset=="Xu")
df_m2 = dplyr::filter(df_all, model=="multiome_powerlaw_v2_noTPMfilter", dataset=="Xu")

if (FALSE) {
	plot_performance_heatmaps(df_m, colors, output_dir, "multiome_powerlaw_v2", comparison=FALSE)
	plot_performance_heatmaps(df_s, colors, output_dir, "scATAC_powerlaw_v2", comparison=FALSE)
	plot_performance_heatmaps(df_m2, colors, output_dir, "multiome_powerlaw_v2_noTPMfilter", comparison=FALSE)
}

## heatmaps of multiome versus scATAC performance
df_comp = make_comparison_df(df_all, stats, model1="multiome_powerlaw_v2", model2="scATAC_powerlaw_v2")
df_comp2  = make_comparison_df(df_all, stats, model1="multiome_powerlaw_v2_noTPMfilter", model2="scATAC_powerlaw_v2")
df_comp3  = make_comparison_df(df_all, stats, model1="multiome_powerlaw_v2", model2="multiome_powerlaw_v2_noTPMfilter")

if (FALSE) {
	plot_performance_heatmaps(df_comp, colors, output_dir, "multiome_minus_scATAC", comparison=TRUE)
	plot_performance_heatmaps(df_comp2, colors, output_dir, "multiome_noTPM_minus_scATAC", comparison=TRUE)
	plot_performance_heatmaps(df_comp3, colors, output_dir, "multiome_minus_multiome_noTPM", comparison=TRUE)

}

# element overlap heatmaps
if (FALSE) {
	plot_jaccard_similarity(jacc, colors, output_dir)
}