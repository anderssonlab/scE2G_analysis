library(tidyverse)
library(cowplot)

proj_dir <- "/oak/stanford/groups/engreitz/Users/sheth/CRISPR_comparison_v3/CRISPR_comparison"
source(file.path(proj_dir, "workflow/scripts/crisprComparisonLoadInputData.R"))
source(file.path(proj_dir, "workflow/scripts/scE2G_bootstrap/crisprComparisonPlotFunctions_trapAUC.R"))
source(file.path(proj_dir, "workflow/scripts/scE2G_bootstrap/crisprComparisonBootstrapFunctions_weighted_trapAUC.R"))

compare_performance_between_categories <- function(merged_training, category_name, category_column, pred_config, models, out_prefix,
    plot_only = FALSE, n_positives_total = 461) {
    ## process inputs pred config
    merged_training <- filter(merged_training, pred_uid %in% models)

    # get colors for all models to plot
    pred_colors <- pred_config %>% 
        filter(pred_uid %in% models) %>% 
        mutate(pred_uid = factor(pred_uid), levels = models, ordered = TRUE) %>% 
        arrange(pred_uid) %>% 
        select(pred_name_long, color) %>% 
        deframe()

    # get thresholds for all predictors
    thresholds <- getThresholdValues(pred_config, threshold_col = "alpha")

    # split input data based on category
    merged_split <- split(merged_training, merged_training[[category_column]])
    n_categories <- length(merged_split)

    # format data for bootstrapping
    merged_bs_split <- lapply(merged_split, convertMergedForBootstrap, pred_config = pred_config)

    ## calculate percent pairs...
    group_counts <- merged_training %>% 
        select(chrom, chromStart, chromEnd, measuredGeneSymbol, Regulated, ExperimentCellType, Dataset, !!sym(category_column)) %>% 
        distinct() %>% mutate(n_total = n()) %>% 
        group_by(!!sym(category_column), n_total) %>% 
        summarize(n_category = n(), n_regulated_in_category = sum(Regulated), .groups = "drop") %>% 
        mutate(n_regulated_total = n_positives_total,
            pct_category_is_regulated = n_regulated_in_category / n_category * 100,
            pct_regulated_in_category = n_regulated_in_category / n_regulated_total * 100,
            pct_pairs_in_category = n_category / n_total * 100,
            category = !!sym(category_column)) %>% 
        arrange(pct_regulated_in_category)

    # category levels
    category_levels <- unique(group_counts$category)
            
    # save group counts
    write_tsv(group_counts, paste0(out_prefix, "_counts.tsv"))

    # reformat group counts
    group_counts <- group_counts %>%
        #pivot_longer(c(pct_regulated_in_category, pct_pairs_in_category), names_to = "metric", values_to = "pct_pairs") %>% 
        #mutate(plot_label = ifelse(metric == "pct_pairs_in_category", "% pairs in category", "% regulated pairs in category")) %>% 
        mutate(category = factor(category, levels = category_levels, ordered = TRUE))
    
    if (plot_only) { message("Plotting only (assuming results exist!)") } else {

        ## calculate auprc and precision
        message("Calculating AUPRC for each category...")
        auprc_res_list <- lapply(merged_bs_split, bootstrapPerformanceIntervals, metric = "auprc", R = 1000)

        #message("Calculating precision for each category...")
        #prec_res_list <- lapply(merged_bs_split, bootstrapPerformanceIntervals, metric = "precision", thresholds = thresholds, R = 1000)

        auprc_res <- rbindlist(auprc_res_list, idcol = "category") %>%
            mutate(metric = "AUPRC") %>%
            select(pred_uid = id, category, metric, full, lower, upper)

        # prec_res <- rbindlist(prec_res_list, idcol = "category") %>%
        #     mutate(metric = "precision") %>%
        #     select(pred_uid = id, category, metric, full, lower, upper)

        if (n_categories == 2) {
            message("Calculating delta AUPRC between categories...")
            delta_auprc <- bootstrapDeltaPerformanceDatasets(merged_bs_split[[1]], merged_bs_split[[2]], metric = "auprc", R = 1000)
            write_tsv(delta_auprc, paste0(out_prefix, "delta_auprc_significance.tsv"))

            # message("Calculating delta precision between categories...")
            # delta_prec <- bootstrapDeltaPerformanceDatasets(merged_bs_split[[1]], merged_bs_split[[2]], metric = "precision", thresholds = thresholds, R = 1000)
            # write_tsv(delta_prec, paste0(out_prefix, "delta_precision_significance.tsv"))
        }

        # combine tables and save
        res <- auprc_res %>% #rbind(auprc_res, prec_res) %>% 
            mutate(pred_name_long = pred_key[pred_uid])

        write_tsv(res, paste0(out_prefix, "performance_summary.tsv"))
        message("Saved performance summary table...")
    }
    
    res <- read_tsv(paste0(out_prefix, "performance_summary.tsv"), show_col_types = FALSE) %>% 
        mutate(pred_name_long = pred_key[pred_uid]) %>% 
        mutate(pred_name_long = factor(pred_name_long, levels = names(pred_colors), ordered = TRUE)) %>%
        mutate(category = factor(category, levels = category_levels, ordered = TRUE))


    ## plot
    # aesthetics
    alpha_key <- setNames((seq(0, 1, length.out = n_categories)), category_levels)
    all_colors <- c(pred_colors,
        pct_pairs_in_category ="#6e788d",
        pct_category_is_regulated = "#435369",
        pct_regulated_in_category = "#435369")
    pct_labs <- c(pct_pairs_in_category = "% tested pairs in category",
        pct_category_is_regulated = "% category that is positive",
        pct_regulated_in_category = "% of positive pairs in category")
    pct_titles <- c(pct_pairs_in_category = "% tested ",
        pct_category_is_regulated = "Positive rate",
        pct_regulated_in_category = "% positives")
    pct_maxes <- c(pct_pairs_in_category = 90, pct_category_is_regulated = 25, pct_regulated_in_category = 100)

    base_theme <- theme_classic() + theme(
            axis.text.x = element_blank(), axis.title.x = element_blank(), 
            axis.text.y = element_text(size = 7, color = "#000000"), axis.title.y = element_text(size = 8),
            plot.title = element_text(color = "#000000", size=10),
            axis.ticks = element_line(color = "#000000"), legend.position = "none")

    pct_list <- list()
    for (i in seq_along(pct_labs)) {
        this_metric <- names(pct_labs)[i]

        pct_list[[this_metric]] <- ggplot(group_counts, aes(x = category, y = !!sym(this_metric))) +
            geom_col(aes(alpha = category), color = all_colors[this_metric], fill = all_colors[this_metric]) +
            scale_alpha_manual(values = alpha_key) + 
            ylim(c(0, pct_maxes[this_metric])) +
            labs(y = pct_labs[this_metric], title = pct_titles[this_metric]) + 
            base_theme
    }
    
    res <- mutate(res, metric_name = ifelse(metric == "AUPRC", "AUPRC", "Precision at threshold"))

    # annotate significance stars per model if delta files exist
    sig_annotations <- NULL
    if (n_categories == 2) {
        delta_files <- c(auprc = paste0(out_prefix, "delta_auprc_significance.tsv"))
        if (all(file.exists(delta_files))) {
            sig_df <- bind_rows(
                read_tsv(delta_files["auprc"], show_col_types = FALSE) %>% mutate(metric = "AUPRC")
            ) %>%
                mutate(pval_label = case_when(
                    pvalue < 0.001 ~ "***",
                    pvalue < 0.01 ~ "**",
                    pvalue < 0.05 ~ "*",
                    TRUE ~ "n.s."
                )) %>%
                filter(pval_label != "") %>%
                mutate(metric_name = ifelse(metric == "AUPRC", "AUPRC", "Precision at threshold"))

            # Keep only one row per model/metric (drop category info)
            sig_annotations <- res %>%
                distinct(pred_uid, pred_name_long, metric, metric_name) %>%
                inner_join(sig_df, by = c("pred_uid" = "id", "metric", "metric_name")) %>%
                mutate(x_stars = 1.5,
                    y_stars = ifelse(metric == "AUPRC", 0.9, 0.8))  # Fixed position near right edge
        }
    }
    
    models <- unique(res$pred_name_long)
    perf_list <- list()
    for (i in seq_along(models)) {
        this_pred <- models[i]
        res_filt <- res %>% filter(pred_name_long == this_pred)

        perf_list[[this_pred]] <- ggplot(res_filt, aes(x = category, y = full)) +
            geom_col(aes(alpha = category), color = all_colors[this_pred], fill = all_colors[this_pred]) +
            geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4, color = "black", alpha = 1) +
            { if (!is.null(sig_annotations)) geom_text(data = filter(sig_annotations, pred_name_long == this_pred),
                aes(x = x_stars, y = y_stars, label = pval_label),
                inherit.aes = FALSE, size = 3.2, hjust = 0.5, vjust = 0.5) else NULL } +
            scale_alpha_manual(values = alpha_key) + 
            ylim(c(0, 1)) +
            labs(y = "AUPRC", title = this_pred) + 
            base_theme
    }

    # save legend
    perf_with_legend <- perf_list[[1]] + theme(legend.position = "right", legend.title = element_text(size = 9), legend.text = element_text(size = 7)) +
        guides(fill = "none", alpha = guide_legend(order = 1, reverse = TRUE))
    legend <- cowplot::get_legend(perf_with_legend)
    ggsave2(paste0(out_prefix, "legend.pdf"), legend, height = 3.5, width = 2)

    # save grid
    if (category_column == "effect_size_group") {
        pct_list[[1]] <- pct_list[[2]]
        pct_list[[2]] <- pct_list[[3]]
    }

    grid <- plot_grid(pct_list[[1]], pct_list[[2]], perf_list[[2]], perf_list[[1]],
        nrow = 1, align = "hv")

    ggsave2(paste0(out_prefix, "performance_summary.pdf"), grid, height = 3, width = 6.5)
    message("Saved performance summary plots!")
}

### MAIN

## other file paths
merged_heldout_crispr_file <- "/oak/stanford/groups/engreitz/Projects/scE2G/data/CRISPR_benchmarking/held_out_250602/expt_pred_merged_annot.txt.gz"
merged_training_crispr_file <- "/oak/stanford/groups/engreitz/Projects/scE2G/data/CRISPR_benchmarking/sc.250425_2.K562.Xu.cv/expt_pred_merged_annot.txt.gz"
pred_config_file <- file.path(proj_dir, "resources", "pred_config", "pred_config_scE2G_heldout.tsv")
correct_annotations_file <- "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_main/2025_0227_CTCF_and_H3K27ac/results/2025_0620_revisit_categories/recategorized_v2/crispr_for_distal_regulation.recategorized.tsv.gz"

cage_annot_file <- "/oak/stanford/groups/engreitz/Projects/scE2G/data/CRISPR_benchmarking/EPCrisprBenchmark_ensemble_data_GRCh38.CAGE_N_0.75.tsv.gz" # CAGE_N_0.75

## load input files
cols_remove <- c("element_category_with_dnase", "element_category_simple",	"element_category", "distance", "distanceToTSS", "Dataset", "data_category")
merged_heldout <- read_tsv(merged_heldout_crispr_file, show_col_types = FALSE) %>%
    select(-any_of(cols_remove))

merged_training <- read_tsv(merged_training_crispr_file, show_col_types = FALSE) %>% 
    select(-any_of(cols_remove)) %>%
    mutate(pred_uid = ifelse(pred_uid == "scE2G_ATAC.E2G.Score.cv.qnorm", "scE2G_ATAC.E2G.Score.qnorm", pred_uid),
        pred_uid = ifelse(pred_uid == "scE2G_multiome.E2G.Score.cv.qnorm", "scE2G_multiome.E2G.Score.qnorm", pred_uid),
        pred_col = ifelse(pred_col == "E2G.Score.cv.qnorm", "E2G.Score.qnorm", pred_col))

pred_config <- importPredConfig(pred_config_file, expr = FALSE)
pred_key <- setNames(pred_config$pred_name_long, pred_config$pred_uid)

## read in correct annotations
cage_annot <- read_tsv(cage_annot_file, show_col_types = FALSE) %>% 
    select(chrom, chromStart, chromEnd, measuredGeneSymbol, Dataset, CAGE_N_0.75)

correct_annot <- read_tsv(correct_annotations_file, show_col_types = FALSE) %>% 
    select(chrom = chr, chromStart = start, chromEnd = end, measuredGeneSymbol, ExperimentCellType = cell_type,
        element_category, direct_vs_indirect_negative, distanceToTSS, Dataset, data_category, ubiq_category)

#correct_annot_heldout <- correct_annot %>% filter(data_category == "validation") %>% left_join(cage_annot, by = c("chrom", "chromStart", "chromEnd", "measuredGeneSymbol", "Dataset"))
correct_annot_training <- correct_annot %>% filter(data_category == "training") %>% left_join(cage_annot, by = c("chrom", "chromStart", "chromEnd", "measuredGeneSymbol", "Dataset"))

## merge with benchmarking data
#merged_heldout <- left_join(merged_heldout, correct_annot_heldout, by = c("chrom", "chromStart", "chromEnd", "measuredGeneSymbol", "ExperimentCellType"))
merged_training <- left_join(merged_training, correct_annot_training, by = c("chrom", "chromStart", "chromEnd", "measuredGeneSymbol", "ExperimentCellType"))

# process merged data for benchmarking analyses, including filtering for ValidConnection == TRUE
merged_training <- processMergedData(merged_training, pred_config = pred_config,
                                     filter_valid_connections = TRUE,
                                     include_missing_predictions = TRUE)

# merged_heldout <- processMergedData(merged_heldout, pred_config = pred_config,
#                                     filter_valid_connections = TRUE,
#                                     include_missing_predictions = TRUE)


### RUN COMPARISONS
out_dir <- file.path(proj_dir, "workflow", "results", "v4_scE2G_category_analysis"); dir.create(out_dir, showWarnings = FALSE)
models <- c("scE2G_multiome.E2G.Score.qnorm", "baseline.distToTSS")

merged_training <- merged_training %>% 
    mutate(chromatin_category = ifelse(element_category == "High H3K27ac", "High H3K27ac", "Other"),
        enhancerness = ifelse(element_category %in% c("High H3K27ac", "H3K27ac"), "H3K27ac+", "Other"),
        distance_category = ifelse(distanceToTSS < 100000, "< 100 kb", ">= 100 kb"),
        CAGE_category = ifelse(CAGE_N_0.75, "High CAGE signal", "Low CAGE signal"))

## promoter class
if (TRUE) {
    category_name <- "Target gene\npromoter class"
    category_column <- "ubiq_category"
    out_prefix <- paste0(out_dir, "/promoter_class_")

    compare_performance_between_categories(merged_training, category_name, category_column, pred_config, models, out_prefix, plot_only = TRUE)
}

## H3K27ac high versus others
if (TRUE) {
    category_name <- "Chromatin state\nof element"
    category_column <- "chromatin_category"
    out_prefix <- paste0(out_dir, "/high_h3k27ac_")

    compare_performance_between_categories(merged_training, category_name, category_column, pred_config, models, out_prefix, plot_only = TRUE)
}

## Enhancernes
if (TRUE) {
    category_name <- "Chromatin state\nof element"
    category_column <- "enhancerness"
    out_prefix <- paste0(out_dir, "/enhancerness_")

    compare_performance_between_categories(merged_training, category_name, category_column, pred_config, models, out_prefix, plot_only = TRUE)
}

## distance
if (TRUE) {
    category_name <- "Distance to TSS"
    category_column <- "distance_category"
    out_prefix <- paste0(out_dir, "/distance_")

    compare_performance_between_categories(merged_training, category_name, category_column, pred_config, models, out_prefix, plot_only = TRUE)
}

## effect size (remove FF option)
if (TRUE) {
    #merged_training_noFF <- merged_training %>% filter(Dataset != "Nasser2021") %>% mutate(effect_size_group = "All (no FlowFISH)"); print(unique(merged_training_noFF$Dataset))
    merged_training_noFF <- merged_training %>% mutate(effect_size_group = "All")
    merged_below5pct <- merged_training_noFF %>% filter(!(Regulated & abs(EffectSize) > 0.05)) %>% mutate(effect_size_group = "Effect size <=5%")
    merged_over5pct <- merged_training_noFF %>% filter(!(Regulated & abs(EffectSize) <= 0.05)) %>% mutate(effect_size_group = "Effect size >5%")
    merged_over10pct <- merged_training_noFF %>% filter(!(Regulated & abs(EffectSize) <= 0.1)) %>% mutate(effect_size_group = "Effect size >10%")

    merged_effect_size <- rbind(merged_below5pct, merged_over5pct)

    category_name <- "Effect size\ngroup"
    category_column <- "effect_size_group"
    out_prefix <- paste0(out_dir, "/effect_size_overunder5")

    compare_performance_between_categories(merged_effect_size, category_name, category_column, pred_config, models, out_prefix, plot_only = TRUE)
}

## Kendall correlation?
if (TRUE) {
    kendall_annot <- merged_training %>% filter(pred_uid == "Kendall.Kendall") %>% select(pair_uid, Kendall_corr = pred_value)
    median_kendall <-  median(kendall_annot$Kendall_corr) # Median Kendall: 0.00568424430753245
    kendall_threshold <- quantile(kendall_annot$Kendall_corr, probs = 0.8) # 80th percentile: 0.02
    message("Kendall threshold: ", kendall_threshold)

    merged_training_kendall <- left_join(merged_training, kendall_annot, by = "pair_uid") %>% 
        mutate(kendall_category = ifelse(Kendall_corr > kendall_threshold, "Above threshold", "Below threshold"))

    category_name <- "Kendall\ncorrelation"
    category_column <- "kendall_category"
    out_prefix <- paste0(out_dir, "/kendall_")

    compare_performance_between_categories(merged_training_kendall, category_name, category_column, pred_config, models, out_prefix, plot_only = TRUE)
}

## CAGE signal at element
if (TRUE) {
    category_name <- "CAGE signal\nat element"
    category_column <- "CAGE_category"
    out_prefix <- paste0(out_dir, "/CAGE_")

    compare_performance_between_categories(merged_training, category_name, category_column, pred_config, models, out_prefix, plot_only = TRUE)
}




