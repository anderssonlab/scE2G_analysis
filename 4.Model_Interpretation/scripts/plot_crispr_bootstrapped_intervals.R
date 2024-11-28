
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



##  MAIN
out_dir <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0710_downsample_for_paper/CRISPR_benchmark_2M"
pred_config <- fread("/oak/stanford/groups/engreitz/Users/sheth/CRISPR_comparison_v2/CRISPR_comparison/resources/pred_config/pred_config_for_bootstrapping.tsv")
model <- "atac" #"multiome" # "atac"

delta_auprc <- fread(file.path(out_dir, paste0(model, "_delta_auprc.tsv")))
delta_precision <- fread(file.path(out_dir, paste0(model, "_delta_precision.tsv")))

df_auprc <- process_bootstrap_data(delta_auprc)
auprc <- plot_intervals(df_auprc, "auprc", pred_config)
ggsave(filename=file.path(out_dir, paste0(model, "_delta_auprc_v1.pdf")), plot=auprc, width=5, height = 4)

df_precision <- process_bootstrap_data(delta_precision)
precision <- plot_intervals(df_precision, "precision", pred_config)
ggsave(filename=file.path(out_dir, paste0(model, "_delta_precision_v1.pdf")), plot=precision, width=5, height = 4)