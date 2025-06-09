suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(ggdist)
  library(ggpubr)
  library(cowplot)
})

is.empty <- function(filename) file.info(filename)$size == 0

get_stats_from_pred_other_models <- function(pred_file, score_column, score_threshold, inverse_predictor, include_promoters) {
	message(pred_file)
    # input files
	pred_full <- fread(pred_file)
    n_all_putative <- nrow(pred_full)

    # threshold predictions
    	# stats from thresholded predictions
    if (inverse_predictor) {
        score_threshold <- -1 * score_threshold
        pred_full <- pred_full %>% mutate(!!sym(score_column) := -1 * !!sym(score_column))
    }

	pred_thresholded <- pred_full %>% dplyr::filter(!!sym(score_column) >= score_threshold)

	# if (!include_promoters) {
	# 	pred_thresholded <- pred_thresholded %>% dplyr::filter(class != "promoter")
	# }

	enh_gene <- pred_thresholded %>% dplyr::select(chr, start, end, TargetGene)

	# get stats
	n_enh <- enh_gene %>% dplyr::select(-c(TargetGene)) %>% distinct() %>% nrow()
	n_gene <- enh_gene %>% dplyr::select(TargetGene) %>% distinct() %>% nrow()
	n_links <- enh_gene %>% nrow()
	
	n_gene_per_enh <- enh_gene %>% group_by(chr, start, end) %>%
		tally() %>% pull(n) %>% mean()
	n_enh_per_gene <- enh_gene %>% group_by(TargetGene) %>% 
		tally() %>% pull(n) %>% mean()

    # distance_cols <- c("distanceToTSS", "distance", "TSS.dist")
    # distance_col_use <- intersect(distance_cols, colnames(pred_thresholded))[1]
    # pred_thresholded$distanceToTSS <- pred_thresholded[[distance_col_use]]

	# mean_distance <- pred_thresholded %>% pull(distanceToTSS) %>% mean()

	mean_size <- enh_gene %>% dplyr::select(chr, start, end) %>% distinct() %>%
		mutate(width = end - start) %>%
		pull(width) %>% mean()

	res <- c(n_enh_elements = n_enh, n_genes_with_enh = n_gene, 
        n_enh_gene_links = n_links, n_putative_links = n_all_putative,
		mean_genes_per_enh = n_gene_per_enh, mean_enh_per_gene = n_enh_per_gene,
		mean_enh_width = mean_size) %>%
		as.list() %>% data.frame()

	return(res)
}

plot_scatter_facet <- function(res, x_var, y_var, cp, label_key, out_dir) {
    if (x_var == "fragments_total") {
        threshold_line <- 2e6
    } else if (x_var == "umi_count") {
        threshold_line <- 1e6
    } else if (x_var == "cell_count") {
        threshold_line <- 100
    } else {threshold_line <- 0}

	# scatter plot
	s <- ggplot(res, aes(x = !!sym(x_var), y = !!sym(y_var), color = pred_name_long)) +
		geom_vline(xintercept=threshold_line, linetype='dashed', color='#c5cad7') + # fragment count ref
		geom_point(size = 2, alpha = 0.75, shape = 16) +
		#stat_cor(aes(label =  after_stat(rr.label)), color = "#000000",  label.x.npc = "left", label.y.npc = "top", size = 3) +
		facet_wrap(vars(pred_name_long), scales = "fixed", axes = "all", ncol = 4) +
		scale_x_log10() + scale_y_log10() +
		scale_color_manual(values = cp, name = "Predictor") +
		labs(x = label_key[x_var], y = label_key[y_var]) + 
		theme_classic() + theme(axis.text = element_text(size = 7, color = "#000000"),
            axis.title = element_text(size = 8), aspect.ratio = 1,
				axis.ticks = element_line(color = "#000000"),
				strip.background = element_blank(), strip.text = element_text(size = 8), legend.position = "none")
    
    out_file <- file.path(out_dir, paste0(y_var, "_by_", x_var, "_noR2.pdf"))
    ggsave(out_file, s, height = 6, width = 9)
}


##### RUN
sample_key <- fread("/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0214_new_global_properties/config/other_model_sample_key.tsv")
pred_config <- fread("/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0214_new_global_properties/config/config_methods.tsv")
include_promoters <- TRUE # pred don't have class column

this_dir <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0214_new_global_properties/prediction_properties"
out_dir <- file.path(this_dir, "suggested_thresholds"); dir.create(out_dir, showWarnings = FALSE)

model_list <- sample_key %>% select(-c(biosample, fragments_total, cell_count, umi_count)) %>% colnames()

# collect properties for each model
if (FALSE) {
    for (model in model_list) {
        message("Working on ", model, "...")
        out_file <- file.path(out_dir, paste0(model, ".tsv"))
        if (file.exists(out_file)) {
            message("Already done.")
        } else {
            this_pred <- sample_key %>% select(biosample, !!sym(model)) %>%
                setNames(c("biosample", "pred_file")) %>%
                filter(!is.na(pred_file), !is.empty(pred_file))
            pred_key <- setNames(this_pred$pred_file, this_pred$biosample)
            this_config <- pred_config %>% filter(method == model)
            
            score_column <- this_config$score_col[1]
            score_threshold <- this_config$threshold[1]
            inverse_predictor <- this_config$inverse_predictor[1]
            out_file <- file.path(out_dir, paste0(model, ".tsv"))
            
            res <- lapply(pred_key, get_stats_from_pred_other_models, score_column, score_threshold, inverse_predictor, include_promoters) %>% 
                rbindlist(idcol = "biosample") %>% 
                as.data.frame() %>% 
                mutate(model_name = model)

            fwrite(res, out_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
        }
    }
}

# combine all models
if (FALSE) {
    out_files <- file.path(out_dir, paste0(model_list, ".tsv"))
    sample_stats <- sample_key %>% select(biosample, fragments_total, cell_count, umi_count)

    res_all <- lapply(out_files, fread) %>% 
        rbindlist() %>% as.data.frame %>% 
        left_join(sample_stats, by = "biosample") %>% 
        left_join(pred_config, by = c("model_name" = "method"))

    fwrite(res_all, file.path(out_dir, "all_models_prediction_properties.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}

# plot
label_key <- c(
    fragments_total = "# unique ATAC fragments in cluster", cell_count = "# cells in cluster", umi_count = "# RNA UMIs in cluster",
    frag_per_cell = "Mean unique ATAC fragments per cell", umi_per_cell = "Mean RNA UMIs per cell",
    n_enh_gene_links = "# enhancer-gene links", n_enh_elements = "# unique enhancers",
    mean_genes_per_enh = "Mean # genes per enhancer", mean_enh_per_gene = "Mean # enhancers per gene",
    n_genes_with_enh = "# genes with 1+ enhancer", n_putative_links = "# candidate enhancer-gene links",
    mean_enh_width = "Mean width of enhancer element (bp)")

if (FALSE) {
    res_all <- fread(file.path(out_dir, "all_models_prediction_properties.tsv")) %>% 
        mutate(pred_name_long = factor(pred_name_long, levels = pred_config$pred_name_long, ordered = TRUE))
    cp <- setNames(pred_config$color, pred_config$pred_name_long)

    for (i in 6:length(label_key)) {
        this_y <- names(label_key)[i]

        plot_scatter_facet(res_all, "fragments_total", this_y, cp, label_key, out_dir)
    }


}

# get summary per model for a property
if (TRUE) {
    this_property <- "n_enh_gene_links"

    res_all <- fread(file.path(out_dir, "all_models_prediction_properties.tsv"))
    smry <- res_all %>% 
        filter(fragments_total > 2e6) %>% 
        select(model_name, pred_name_long, prop = !!sym(this_property)) %>% 
        filter(prop > 0) %>%
        group_by(model_name, pred_name_long) %>% 
        summarize(min = min(prop), max = max(prop), mean = mean(prop), median = median(prop)) %>% 
        mutate(oom_var = abs(log10(max) - log10(min)))

    fwrite(smry, file.path(out_dir, paste0("summary_of_", this_property, ".tsv")), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

}