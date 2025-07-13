
# libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(data.table)
  library(cowplot)
  library(ggdist)
  library(stringr)
})



plot_crispr_performance <- function(perf, out_dir) {
	# column names: cell_type	pred_uid	pred_name_long	AUPRC	AUPRC_lowerCi	AUPRC_upperCi	PrecThresh	PrecThresh_lowerCi	PrecThresh_upperCi	threshold	PrecMinSens	PrecMinSens_lowerCi	PrecMinSens_upperCi	thresholdMinSens

	# get order of predictors (based on auprc)
	scE2G <- c("scE2G (Multiome)", "scE2G (Multiome), 100 cells x 20k fragments", "scE2G (Multiome), 200 cells x 10k fragments", "scE2G (Multiome), 500 cells x 4k fragments",
		"scE2G (ATAC)", "scE2G (ATAC), 100 cells x 20k fragments", "scE2G (ATAC), 200 cells x 10k fragments", "scE2G (ATAC), 500 cells x 4k fragments")

	order_df <- dplyr::select(perf, pred_uid, pred_name_long, AUPRC) %>%
		dplyr::filter(!str_detect(pred_name_long, "^scE2G")) %>%
		arrange(desc(AUPRC))
	pred_order <- c(scE2G, order_df$pred_name_long)

	# make "grouping" variable for plotting
	perf <- perf %>%
		mutate(perf, 
		pred_name_long = factor(pred_name_long, levels = pred_order, ordered = TRUE),
		grouping = case_when(str_detect(pred_name_long, "scE2G \\(Multiome") ~ "aaa",
			str_detect(pred_name_long, "scE2G \\(ATAC") ~ "bbb",
			.default = "ccc"))

	pred_colors <- setNames(perf$color, perf$pred_name_long)

	auprc <- ggplot(perf, aes(x = pred_name_long, y = AUPRC, fill = pred_name_long)) +
    	geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    	geom_errorbar(aes(ymin = AUPRC_lowerCi, ymax = AUPRC_upperCi), position = position_dodge(width = 0.9),
                  color = "black", width = 0.25) +
		facet_grid(.~grouping, scale="free_x", space="free") +
    	labs(fill = "Predictor", x = NULL, y = "AUPRC") +
    	scale_fill_manual(values = pred_colors[levels(perf$pred_name_long)]) + 
    	scale_y_continuous(limits = c(0, 1)) +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
			axis.ticks.x = element_blank(), axis.text.x  = element_blank(), strip.background = element_blank(),
  			strip.text.x = element_blank(), legend.position = "right")

	prec <- ggplot(perf, aes(x = pred_name_long, y = PrecThresh, fill = pred_name_long)) +
    	geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    	geom_errorbar(aes(ymin = PrecThresh_lowerCi, ymax = PrecThresh_upperCi), position = position_dodge(width = 0.9),
                  color = "black", width = 0.25) +
		facet_grid(.~grouping, scale="free_x", space="free") +
    	labs(fill = "Predictor", x = NULL, y = "Precision at model threshold") +
    	scale_fill_manual(values = pred_colors[levels(perf$pred_name_long)]) + 
    	scale_y_continuous(limits = c(0, 1)) +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
			axis.ticks.x = element_blank(), axis.text.x  = element_blank(), strip.background = element_blank(),
  			strip.text.x = element_blank(), legend.position = "None")

	legend_v <- as_grob(cowplot::get_plot_component(auprc, "guide-box-right"))
	auprc <- auprc + theme(legend.position = "None")

	col1 <- plot_grid(auprc, prec, nrow = 2, ncol = 1, rel_heights = c(1, 1))
	full <- plot_grid(col1, legend_v, nrow = 1, ncol = 2, rel_widths = c(1.3, 1))
	
	out_file <- file.path(out_dir, "CRISPR_performance_all_v1.pdf")
	ggsave2(filename=out_file, plot = full, width= 9, height = 6)
}

##  MAIN
input_dir <- ("/oak/stanford/groups/engreitz/Users/sheth/CRISPR_comparison_v3/CRISPR_comparison/workflow/results/2025_0421_scE2G_downsample_vs_others")
out_dir <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0124_new_downsample/CRISPR_benchmark_2M"; dir.create(out_dir, showWarnings = FALSE)
auprc_file <- file.path(input_dir, "all_pred_auprc_ci.tsv")
prec_file <- file.path(input_dir, "all_pred_precision_ci.tsv")
pred_config <- fread("/oak/stanford/groups/engreitz/Users/sheth/CRISPR_comparison_v3/CRISPR_comparison/resources/pred_config/pred_config_scE2G.tsv")
pred_exclude <- c("ARC.ARC.E2G.Score")

pred_config <- mutate(pred_config, pred_uid = paste0(pred_id, ".", pred_col))

auprc <- fread(auprc_file) %>% select(id, AUPRC = full, AUPRC_lowerCi = lower, AUPRC_upperCi = upper)
prec <- fread(prec_file) %>% select(id, PrecThresh = full, PrecThresh_lowerCi = lower, PrecThresh_upperCi = upper)

perf <- left_join(auprc, prec, by = "id") %>%
	rename(pred_uid = id) %>% 
	filter(!(pred_uid %in% pred_exclude)) %>%
	left_join(pred_config, by = "pred_uid")

a <- plot_crispr_performance(perf, out_dir)
