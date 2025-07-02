library(tidyverse)
library(cowplot)

proj_dir <- "/oak/stanford/groups/engreitz/Users/sheth/CRISPR_comparison_v3/CRISPR_comparison"
source(file.path(proj_dir, "workflow/scripts/crisprComparisonLoadInputData.R"))
source(file.path(proj_dir, "workflow/scripts/scE2G_bootstrap/crisprComparisonPlotFunctions_trapAUC.R"))
source(file.path(proj_dir, "workflow/scripts/scE2G_bootstrap/crisprComparisonBootstrapFunctions_weighted_trapAUC.R"))

compare_performance_between_categories <- function(merged_training, category_name, category_column, pred_config, models, out_prefix, plot_only = FALSE) {
    ## process inputs pred config
    merged_training <- filter(merged_training, pred_uid %in% models)

    # category levels
    category_levels <- unique(merged_training[[category_column]])
    category_levels <- rev(category_levels)


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
        mutate(pct_regulated_in_category = n_regulated_in_category / n_category * 100,
            pct_pairs_in_category = n_category / n_total * 100,
            category = !!sym(category_column)) %>%
        pivot_longer(c(pct_regulated_in_category, pct_pairs_in_category), names_to = "metric", values_to = "pct_pairs") %>% 
        mutate(plot_label = ifelse(metric == "pct_pairs_in_category", "% pairs in category", "% regulated pairs in category")) %>% 
        mutate(category = factor(category, levels = category_levels, ordered = TRUE))
    
    if (plot_only) { message("Plotting only (assuming results exist!)") } else {

        ## calculate auprc and precision
        message("Calculating AUPRC for each category...")
        auprc_res_list <- lapply(merged_bs_split, bootstrapPerformanceIntervals, metric = "auprc", R = 1000)

        message("Calculating precision for each category...")
        prec_res_list <- lapply(merged_bs_split, bootstrapPerformanceIntervals, metric = "precision", thresholds = thresholds, R = 1000)

        auprc_res <- rbindlist(auprc_res_list, idcol = "category") %>%
            mutate(metric = "AUPRC") %>%
            select(pred_uid = id, category, metric, full, lower, upper)

        prec_res <- rbindlist(prec_res_list, idcol = "category") %>%
            mutate(metric = "precision") %>%
            select(pred_uid = id, category, metric, full, lower, upper)

        if (n_categories == 2) {
            message("Calculating delta AUPRC between categories...")
            delta_auprc <- bootstrapDeltaPerformanceDatasets(merged_bs_split[[1]], merged_bs_split[[2]], metric = "auprc", R = 1000)
            write_tsv(delta_auprc, paste0(out_prefix, "delta_auprc_significance.tsv"))

            message("Calculating delta precision between categories...")
            delta_prec <- bootstrapDeltaPerformanceDatasets(merged_bs_split[[1]], merged_bs_split[[2]], metric = "precision", thresholds = thresholds, R = 1000)
            write_tsv(delta_prec, paste0(out_prefix, "delta_precision_significance.tsv"))
        }

        # combine tables and save
        res <- rbind(auprc_res, prec_res) %>% 
            mutate(pred_uid = factor(pred_uid, levels = models, ordered = TRUE),
                pred_name_long = pred_key[pred_uid]) %>% 
            arrange(pred_uid) %>% 
            mutate(pred_name_long = factor(pred_name_long, levels = names(pred_colors), ordered = TRUE))
        write_tsv(res, paste0(out_prefix, "performance_summary.tsv"))
        message("Saved performance summary table...")
    }
    
    res <- read_tsv(paste0(out_prefix, "performance_summary.tsv")) %>% 
        mutate(pred_name_long = pred_key[pred_uid]) %>% 
        mutate(pred_name_long = factor(pred_name_long, levels = names(pred_colors), ordered = TRUE)) %>%
        mutate(category = factor(category, levels = category_levels, ordered = TRUE))


    ## plot
    # aesthetics
    alpha_key <- setNames((seq(0.4, 1, length.out = n_categories)), rev(unique(merged_training[[category_column]])))
    cat_key <- c(`% pairs in category` ="#96a0b3", `% regulated pairs in category` = "#435369")

    pct_pairs <- ggplot(group_counts, aes(x = pct_pairs, y = category)) +
        #geom_col(aes(fill = plot_label)) +
        geom_col(aes(alpha = category), fill = "#435369") +
        facet_wrap(~ plot_label, ncol = 1, strip.position = "bottom", scales = "free") +
        scale_alpha_manual(values = alpha_key) + #scale_fill_manual(values = cat_key) +
        labs(y = "Element-gene pair category") + 
        theme_classic() + theme(strip.placement = "outside", strip.background = element_blank(), axis.title.x = element_blank(), # facet label is x axis title
            axis.text = element_text(size = 7, color = "#000000"), axis.title.y = element_text(size = 9),
            axis.ticks = element_line(color = "#000000"), legend.position = "none")

    res <- mutate(res, metric_name = ifelse(metric == "AUPRC", "AUPRC", "Precision at threshold"))

    # annotate significance stars per model if delta files exist
    sig_annotations <- NULL
    if (n_categories == 2) {
        delta_files <- c(auprc = paste0(out_prefix, "delta_auprc_significance.tsv"),
                        precision = paste0(out_prefix, "delta_precision_significance.tsv"))
        if (all(file.exists(delta_files))) {
            sig_df <- bind_rows(
                read_tsv(delta_files["auprc"], show_col_types = FALSE) %>% mutate(metric = "AUPRC"),
                read_tsv(delta_files["precision"], show_col_types = FALSE) %>% mutate(metric = "precision")
            ) %>%
                mutate(pval_label = case_when(
                    pvalue < 0.001 ~ "***",
                    pvalue < 0.01 ~ "**",
                    pvalue < 0.05 ~ "*",
                    TRUE ~ ""
                )) %>%
                filter(pval_label != "") %>%
                mutate(metric_name = ifelse(metric == "AUPRC", "AUPRC", "Precision at threshold"))

            # Keep only one row per model/metric (drop category info)
            sig_annotations <- res %>%
                distinct(pred_uid, pred_name_long, metric, metric_name) %>%
                inner_join(sig_df, by = c("pred_uid" = "id", "metric", "metric_name")) %>%
                mutate(y = as.numeric(pred_name_long),
                    x = ifelse(metric == "AUPRC", 0.85, 0.8))  # Fixed position near right edge
        }
    }

    perf <- ggplot(res, aes(y = pred_name_long, x = full, fill = pred_name_long, alpha = category)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_errorbar(aes(xmin = lower, xmax = upper, group = category),
                        position = position_dodge(0.9), width = 0.4, color = "black", alpha = 1) +
        { if (!is.null(sig_annotations)) geom_text(data = sig_annotations,
            aes(x = x, y = y, label = pval_label),
            inherit.aes = FALSE, size = 3.2, hjust = 1, vjust = 0.5) else NULL } +
        facet_wrap(~ metric_name, ncol = 2, strip.position = "bottom") +
        labs(alpha = category_name) +
        scale_fill_manual(values = pred_colors, na.value = NA) + scale_alpha_manual(values = alpha_key) +
        scale_x_continuous(limits = c(0, 1)) +
        guides(fill = "none") +
        theme_classic() + theme(strip.placement = "outside", strip.background = element_blank(), axis.title = element_blank(), # facet label is x axis title
            axis.text = element_text(size = 7, color = "#000000"), 
            axis.ticks = element_line(color = "#000000"), legend.position = "none")

    perf_with_legend <- perf + theme(legend.position = "right", legend.title = element_text(size = 9), legend.text = element_text(size = 7)) +
        guides(fill = "none", alpha = guide_legend(order = 1, reverse = TRUE))
    legend <- cowplot::get_legend(perf_with_legend)

    left <- plot_grid(NULL, pct_pairs, NULL, ncol = 1, rel_heights = c(0.75, 3, 0.75))
    grid <- plot_grid(pct_pairs, perf, legend, nrow = 1, rel_widths = c(0.95, 1.55, 0.4))
    ggsave2(paste0(out_prefix, "performance_summary.pdf"), height = 3.5 * 2/3, width = 9)
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

pred_config <- importPredConfig(pred_config_file, expr = FALSE)
pred_key <- setNames(pred_config$pred_name_long, pred_config$pred_uid)

## read in correct annotations
cage_annot <- fread(cage_annot_file) %>% 
    select(chrom, chromStart, chromEnd, measuredGeneSymbol, Dataset, CAGE_N_0.75)

correct_annot <- read_tsv(correct_annotations_file) %>% 
    select(chrom = chr, chromStart = start, chromEnd = end, measuredGeneSymbol, ExperimentCellType = cell_type,
        element_category, direct_vs_indirect_negative, distanceToTSS, Dataset, data_category)

correct_annot_heldout <- correct_annot %>% filter(data_category == "validation") %>% 
    left_join(cage_annot, by = c("chrom", "chromStart", "chromEnd", "measuredGeneSymbol", "Dataset"))
correct_annot_training <- correct_annot %>% filter(data_category == "training") %>% 
    left_join(cage_annot, by = c("chrom", "chromStart", "chromEnd", "measuredGeneSymbol", "Dataset"))

## merge with benchmarking data
merged_heldout <- left_join(merged_heldout, correct_annot_heldout, by = c("chrom", "chromStart", "chromEnd", "measuredGeneSymbol", "ExperimentCellType"))
merged_training <- left_join(merged_training, correct_annot_training, by = c("chrom", "chromStart", "chromEnd", "measuredGeneSymbol", "ExperimentCellType"))

# process merged data for benchmarking analyses, including filtering for ValidConnection == TRUE
merged_training <- processMergedData(merged_training, pred_config = pred_config,
                                     filter_valid_connections = TRUE,
                                     include_missing_predictions = TRUE)

# merged_heldout <- processMergedData(merged_heldout, pred_config = pred_config,
#                                     filter_valid_connections = TRUE,
#                                     include_missing_predictions = TRUE)


### RUN COMPARISONS
out_dir <- file.path(proj_dir, "workflow", "results", "v3_scE2G_category_analysis"); dir.create(out_dir, showWarnings = FALSE)
models <- c("scE2G_multiome.E2G.Score.qnorm", "scE2G_ATAC.E2G.Score.qnorm", "scATAC_ABC.ABC.Score", "ARC.ARC.E2G.Score", "Kendall.Kendall")

merged_training <- merged_training %>% 
    mutate(chromatin_category = ifelse(element_category == "High H3K27ac", "High H3K27ac", "Other"),
        enhancerness = ifelse(element_category %in% c("High H3K27ac", "H3K27ac"), "H3K27ac+", "Other"),
        distance_category = ifelse(distanceToTSS < 100000, "< 100 kb", ">= 100 kb"),
        CAGE_category = ifelse(CAGE_N_0.75, "High CAGE signal", "Low CAGE signal"))

## promoter class
if (FALSE) {
    category_name <- "Target gene\npromoter class"
    category_column <- "ubiq_category"
    out_prefix <- paste0(out_dir, "/promoter_class_")

    compare_performance_between_categories(merged_training, category_name, category_column, pred_config, models, out_prefix, plot_only = FALSE)
}

## H3K27ac high versus others
if (TRUE) {
    category_name <- "Chromatin state\nof element"
    category_column <- "chromatin_category"
    out_prefix <- paste0(out_dir, "/high_h3k27ac_")

    compare_performance_between_categories(merged_training, category_name, category_column, pred_config, models, out_prefix, plot_only = FALSE)
}

## Enhancernes
if (TRUE) {
    category_name <- "Chromatin state\nof element"
    category_column <- "enhancerness"
    out_prefix <- paste0(out_dir, "/enhancerness_")

    compare_performance_between_categories(merged_training, category_name, category_column, pred_config, models, out_prefix, plot_only = FALSE)
}

## distance
if (FALSE) {
    category_name <- "Distance to TSS"
    category_column <- "distance_category"
    out_prefix <- paste0(out_dir, "/distance_")

    compare_performance_between_categories(merged_training, category_name, category_column, pred_config, models, out_prefix, plot_only = FALSE)
}

## effect size (remove FF option)
if (FALSE) {
    merged_training_noFF <- merged_training %>% filter(Dataset != "Nasser2021") %>% mutate(effect_size_group = "All (no FlowFISH)"); print(unique(merged_training_noFF$Dataset))
    merged_over5pct <- merged_training_noFF %>% filter(!(Regulated & abs(EffectSize) <= 0.05)) %>% mutate(effect_size_group = "Effect size >5%")
    merged_over10pct <- merged_training_noFF %>% filter(!(Regulated & abs(EffectSize) <= 0.1)) %>% mutate(effect_size_group = "Effect size >10%")

    merged_effect_size <- rbind(merged_training_noFF, merged_over5pct, merged_over10pct)

    category_name <- "Effect size\ngroup"
    category_column <- "effect_size_group"
    out_prefix <- paste0(out_dir, "/effect_size_")

    compare_performance_between_categories(merged_effect_size, category_name, category_column, pred_config, models, out_prefix)
}

## Kendall correlation?
if (FALSE) {
    kendall_annot <- merged_training %>% filter(pred_uid == "Kendall.Kendall") %>% select(pair_uid, Kendall_corr = pred_value)
    median_kendall <-  median(kendall_annot$Kendall_corr) # Median Kendall: 0.00568424430753245
    kendall_threshold <- quantile(kendall_annot$Kendall_corr, probs = 0.8) # 80th percentile: 0.02
    message("Kendall threshold: ", kendall_threshold)

    merged_training_kendall <- left_join(merged_training, kendall_annot, by = "pair_uid") %>% 
        mutate(kendall_category = ifelse(Kendall_corr > kendall_threshold, "Above threshold", "Below threshold"))

    category_name <- "Kendall\ncorrelation"
    category_column <- "kendall_category"
    out_prefix <- paste0(out_dir, "/kendall_")

    compare_performance_between_categories(merged_training_kendall, category_name, category_column, pred_config, models, out_prefix, plot_only = FALSE)
}

## CAGE signal at element
if (FALSE) {
    category_name <- "CAGE signal\nat element"
    category_column <- "CAGE_category"
    out_prefix <- paste0(out_dir, "/CAGE_")

    compare_performance_between_categories(merged_training, category_name, category_column, pred_config, models, out_prefix, plot_only = FALSE)
}




