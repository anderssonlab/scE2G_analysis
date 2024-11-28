
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

process_bootstrap_data <- function(df) {
	# df header: id	metric	full	conf	lower	upper	min	max	pvalue

	df <- df %>% separate_wider_delim(id, delim = " | ", names = c("model1", "model2")) %>%
		dplyr::filter(str_detect(model1, "^scE2G") + str_detect(model2, "^scE2G") == 1)  %>% # select rows with exactly one scE2G model
		mutate(sign_label = case_when(pvalue < 0.001 ~ "***",
			pvalue < 0.01 ~ "**",
			pvalue < 0.05 ~ "*",
			pvalue >= 0.05 ~ " "))

	e2g_first <- dplyr::filter(df, str_detect(model1, "^scE2G")) %>%
		mutate(model_temp = model2, model2 = model1, model1 = model_temp,
		 full = -full, old_lower = lower, lower = -upper, upper = -old_lower,
		 old_min = min, min = -max, max = -old_min) %>%
		dplyr::select(-c(model_temp, old_lower, old_min))
	e2g_second <- dplyr::filter(df, str_detect(model2, "^scE2G"))
		
	df <- rbind(e2g_first, e2g_second) %>%
		mutate(label_pos = ifelse(lower < 0, lower - 0.03, lower + 0.03))

	return(df)
}

plot_intervals <- function(df, metric, pred_config) {
	pred_config <- pred_config %>% mutate(pred_uid = paste0(pred_id, ".", pred_col))
	pred_colors <- setNames(pred_config$color, pred_config$pred_id)

	# get order of predictors
	ord <- df %>% group_by(model1) %>%
		summarize(avg_delta = mean(full)) %>%
		arrange(avg_delta) %>%
		pull(model1)

	df <- df %>% mutate(comp_id = paste0(model1, ".", model2)) %>%
		left_join(select(pred_config, pred_uid, pred_id), by = c("model1"="pred_uid")) %>%
		mutate(model1 = factor(model1, levels = ord, ordered = TRUE)) %>%
		arrange(model1) %>%
		mutate(comp_id = factor(comp_id, levels = comp_id, ordered = TRUE))

	df$pred_id = factor(df$pred_id, levels = unique(df$pred_id), ordered = TRUE)

	# some params
	x_label <- ifelse(metric == "auprc", "Delta AUPRC from scE2G on 2M fragment cluster", "Delta precision from scE2G on 2M fragment cluster")
	x_min <- min(min(df$label_pos), min(df$lower)) - 0.05
	x_max <- max(max(df$label_pos), max(df$upper))

	g <- ggplot(df, aes(x = full, y = comp_id, fill = pred_id)) + 
		geom_bar(stat = "identity", width = 0.8, position  = "dodge") +
		geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.3, linewidth = 0.25) +
		geom_text(aes(label = sign_label, x = label_pos), color = "#000000", size = 2) +
		scale_fill_manual(values = pred_colors, guide = guide_legend(reverse = TRUE)) +
		xlab(x_label) + ylab("Model") +
		#xlim(c(x_min, x_max))+
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), axis.text.y = element_blank(),
  			legend.position="right", legend.direction="vertical", legend.text=element_text(size=7), legend.title=element_blank())
	
	return(g)
}

plot_crispr_performance <- function(perf, out_dir) {
	# column names: cell_type	pred_uid	pred_name_long	AUPRC	AUPRC_lowerCi	AUPRC_upperCi	PrecThresh	PrecThresh_lowerCi	PrecThresh_upperCi	threshold	PrecMinSens	PrecMinSens_lowerCi	PrecMinSens_upperCi	thresholdMinSens

	# get order of predictors (based on auprc)
	scE2G <- c("scE2G (Multiome)", "scE2G (Multiome), 2k cells x 1k fragments", "scE2G (Multiome), 1k cells x 2k fragments", "scE2G (Multiome), 500 cells x 4k fragments",
		"scE2G (ATAC)", "scE2G (ATAC), 2k cells x 1k fragments", "scE2G (ATAC), 1k cells x 2k fragments", "scE2G (ATAC), 500 cells x 4k fragments")
	# c("scE2G (Multiome,2000_cells_1000_atac)", "scE2G (Multiome,1000_cells_2000_atac)", "scE2G (Multiome,500_cells_4000_atac)",
	# 	"scE2G (ATAC,2000_cells_1000_atac)", "scE2G (ATAC,1000_cells_2000_atac)", "scE2G (ATAC,500_cells_4000_atac)")
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
out_dir <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0710_downsample_for_paper/CRISPR_benchmark_2M"
input_file <- file.path(out_dir, "performance_summary.txt")
prec_file <- file.path(out_dir, "all_pred_precision_ci.tsv")
pred_config <- fread("/oak/stanford/groups/engreitz/Users/sheth/CRISPR_comparison_v2/CRISPR_comparison/resources/pred_config/pred_config_for_bootstrapping.tsv")
model <- "atac" #"multiome" # "atac"

pred_config <- mutate(pred_config, pred_uid = paste0(pred_id, ".", pred_col))

prec <- fread(prec_file) %>%
	setNames(c("pred_uid", "metric", "PrecThresh", "conf", "PrecThresh_lowerCi", "PrecThresh_upperCi", "PrecThresh_min", "PrecThresh_max")) %>%
	dplyr::select(pred_uid, PrecThresh, PrecThresh_lowerCi, PrecThresh_upperCi)

perf <- fread(input_file, sep = "\t") %>%
	dplyr::filter(pred_uid %in% pred_config$pred_uid) %>%
	dplyr::select(-c(pred_name_long, PrecThresh, PrecThresh_lowerCi, PrecThresh_upperCi)) %>%
	left_join(pred_config, by = "pred_uid") %>%
	left_join(prec, by = "pred_uid")

a <- plot_crispr_performance(perf, out_dir)
