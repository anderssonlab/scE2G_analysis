# libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(data.table)
  library(cowplot)
  library(ggdist)
  library(scales)
})

plot_crispr_effect_size <-function(results_dir, dataset, model_name, threshold, outdir){
	score_name = "E2G.Score.cv.qnorm"
	crispr_file = file.path(results_dir, dataset, model_name, "EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz")
	outfile = file.path(outdir, paste0(dataset, ".", model_name, ".pdf"))
	crispr_data = fread(crispr_file, sep="\t", header=TRUE)

	crispr_data$model_score = crispr_data[[score_name]]
	crispr_data$EffectSize = as.numeric(crispr_data$EffectSize) * 100

	# add category
	crispr_data$category = "empty"
	crispr_data$category[crispr_data$Significant==FALSE] = "Not significant"
	crispr_data$category[crispr_data$Significant==TRUE] = "Increase"
	crispr_data$category[crispr_data$Regulated==TRUE] = "Decrease"
	

	crispr_data = dplyr::filter(crispr_data, model_score>0) %>%
		arrange(desc(category)) # not sig, increase, decrease 

	cp = c("#6e788d", "#c5373d", "#006eae")
	names(cp) = c("Not significant", "Decrease", "Increase")

	g = ggplot(crispr_data, aes(x=model_score, y=EffectSize, color=category)) +
		geom_point(alpha=0.7) +
		geom_vline(xintercept=threshold, linetype="dashed", color="#1c2a43") +
		xlab("scE2G score") + ylab("Perturbation effect size (%)") +
		scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
		scale_color_manual(values=cp) +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position="right", legend.title=element_blank()) 

	ggsave(file=outfile, plot=g, width=5, height=4)
}

results_dir = "/oak/stanford/groups/engreitz/Users/sheth/sc-E2G/results/2024_0826_K562_GM12878"
output_dir = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0710_downsample_for_paper/crispri_effect_size"

dir.create(output_dir)

plot_crispr_effect_size(results_dir, dataset="K562_Xu_pl", model_name="multiome_powerlaw_v2", threshold=0.164, outdir = output_dir)
plot_crispr_effect_size(results_dir, dataset="K562_Xu_pl", model_name="scATAC_powerlaw_v2", threshold=0.183, outdir = output_dir)


