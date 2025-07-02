library(tidyverse)
library(cowplot)

proj_dir <- "/oak/stanford/groups/engreitz/Users/sheth/CRISPR_comparison_v3/CRISPR_comparison"
source(file.path(proj_dir, "workflow/scripts/crisprComparisonLoadInputData.R"))
source(file.path(proj_dir, "workflow/scripts/scE2G_bootstrap/crisprComparisonPlotFunctions_trapAUC.R"))
source(file.path(proj_dir, "workflow/scripts/scE2G_bootstrap/crisprComparisonBootstrapFunctions_weighted_trapAUC.R"))

## file paths
merged_heldout_crispr_file <- "/oak/stanford/groups/engreitz/Projects/scE2G/data/CRISPR_benchmarking/held_out_250602/expt_pred_merged_annot.txt.gz"
merged_training_crispr_file <- "/oak/stanford/groups/engreitz/Projects/scE2G/data/CRISPR_benchmarking/sc.250425_2.K562.Xu.cv/expt_pred_merged_annot.txt.gz"
pred_config_file <- file.path(proj_dir, "resources", "pred_config", "pred_config_scE2G_heldout.tsv")
correct_annotations_file <- "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_main/2025_0227_CTCF_and_H3K27ac/results/2025_0620_revisit_categories/recategorized_v2/crispr_for_distal_regulation.recategorized.tsv.gz"

## load input files
cols_remove <- c("element_category_with_dnase", "element_category_simple",	"element_category", "distance", "distanceToTSS", "Dataset", "data_category")
merged_heldout <- read_tsv(merged_heldout_crispr_file, show_col_types = FALSE) %>%
    select(-any_of(cols_remove), -starts_with("direct_"), -starts_with("indirect_"))

merged_training <- read_tsv(merged_training_crispr_file, show_col_types = FALSE) %>% 
    select(-any_of(cols_remove), -tarts_with("direct_"), -starts_with("indirect_")) %>%
    mutate(pred_uid = ifelse(pred_uid == "scE2G_ATAC.E2G.Score.cv.qnorm", "scE2G_ATAC.E2G.Score.qnorm", pred_uid),
        pred_uid = ifelse(pred_uid == "scE2G_multiome.E2G.Score.cv.qnorm", "scE2G_multiome.E2G.Score.qnorm", pred_uid),
        pred_col = ifelse(pred_col == "E2G.Score.cv.qnorm", "E2G.Score.qnorm", pred_col))

pred_config <- importPredConfig(pred_config_file, expr = FALSE)
pred_key <- setNames(pred_config$pred_name_long, pred_config$pred_uid)

#### SET UP HERE ####
## filter to specific chromatin categoreis and/or effect sizes, distance to tss
out_dir <- file.path(proj_dir, "workflow", "results", "v4_scE2G_heldout_analysis"); dir.create(out_dir, showWarnings = FALSE)
out_dir <- file.path(out_dir, "noK562_filter_CTCF_H3K27me3_noH3K27ac_ES5pct");  dir.create(out_dir, showWarnings = FALSE)
#out_dir <- file.path(out_dir, "no_filter");  dir.create(out_dir, showWarnings = FALSE)
plot_only <- FALSE
w <- 12; h <- 8

correct_annot <- read_tsv(correct_annotations_file) %>% 
    select(chrom = chr, chromStart = start, chromEnd = end, measuredGeneSymbol, ExperimentCellType = cell_type,
        element_category, direct_vs_indirect_negative, distanceToTSS, Dataset, data_category)

correct_annot_heldout <- correct_annot %>% filter(data_category == "validation")
correct_annot_training <- correct_annot %>% filter(data_category == "training")

merged_heldout <- left_join(merged_heldout, correct_annot_heldout, by = c("chrom", "chromStart", "chromEnd", "measuredGeneSymbol", "ExperimentCellType")) %>% 
    filter(distanceToTSS < 1000000) %>% 
    filter(!(Regulated & abs(EffectSize) <= 0.05)) %>% 
    filter(!(Regulated & element_category %in% c("CTCF element", "H3K27me3 element"))) %>% 
    filter(!(Regulated & element_category == "No H3K27ac")) %>% 
    filter(ExperimentCellType != "K562")

merged_training <- left_join(merged_training, correct_annot_training, by = c("chrom", "chromStart", "chromEnd", "measuredGeneSymbol", "ExperimentCellType")) %>% 
    filter(distanceToTSS < 1000000) %>% 
    filter(!(Regulated & abs(EffectSize) <= 0.05)) %>% 
    filter(!(Regulated & element_category %in% c("CTCF element", "H3K27me3 element"))) %>% 
    filter(!(Regulated & element_category == "No H3K27ac"))

## summarize number of positives
heldout_summary <- merged_heldout %>%
    select(chrom, chromStart, chromEnd, measuredGeneSymbol, Regulated, ExperimentCellType, Dataset, data_category, direct_vs_indirect_negative) %>% 
    distinct()
training_summary <- merged_training %>% 
    select(chrom, chromStart, chromEnd, measuredGeneSymbol, Regulated, ExperimentCellType, Dataset, data_category, direct_vs_indirect_negative) %>% 
    distinct()
combined_summary <- rbind(heldout_summary, training_summary) %>% 
    group_by(ExperimentCellType, Dataset, data_category) %>% 
    summarize(n_total = n(), n_positive = sum(Regulated),
        weighted_positives = sum(Regulated * direct_vs_indirect_negative))
write_tsv(combined_summary, file.path(out_dir, "crispr_dataset_summary.tsv"))

category_summary <- combined_summary %>%
    group_by(data_category) %>%
    summarize(n_total = sum(n_total), n_positive = sum(n_positive), weighted_positives = sum(weighted_positives))
write_tsv(category_summary, file.path(out_dir, "crispr_dataset_category_summary.tsv"))

n_key <- category_summary %>% select(data_category, n_positive) %>% deframe()
n_weighted_key <- category_summary %>% mutate(weighted_positives = round(weighted_positives, 2)) %>%
    select(data_category, weighted_positives) %>% deframe()
unweighted_title <- paste0("# positives, training: ", n_key["training"], "\n", "# positives, held-out: ", n_key["validation"])
weighted_title <- paste0("Weighted positives, training: ", n_weighted_key["training"], "\n", "Weighted positives, held-out: ", n_weighted_key["validation"])
weighted_heldout_title <- paste0("Weighted positives, held-out: ", n_weighted_key["validation"])

# process merged data for benchmarking analyses, including filtering for ValidConnection == TRUE
merged_training <- processMergedData(merged_training, pred_config = pred_config,
                                     filter_valid_connections = TRUE,
                                     include_missing_predictions = TRUE)

merged_heldout <- processMergedData(merged_heldout, pred_config = pred_config,
                                    filter_valid_connections = TRUE,
                                    include_missing_predictions = TRUE)

# get models available for both the training and held out CRISPR data
models <- intersect(merged_training$pred_uid, merged_heldout$pred_uid) %>% unique

# only retain data on models to include in analyses
merged_training <- filter(merged_training, pred_uid %in% models)
merged_heldout <- filter(merged_heldout, pred_uid %in% models)

# get colors for all models to plot
pred_colors <- pred_config %>% 
  filter(pred_uid %in% models) %>% 
  select(pred_name_long, color) %>% 
  deframe()

# get fill colors, which includes information on CRISPR dataset type
fill_colors <- c(Training = NA_character_, pred_colors)

# ggplot themes
base_theme <- theme(text = element_text(size = 13, color = "#000000"),
    axis.text = element_text(color = "#000000"),
    axis.ticks = element_line(color = "#000000"))
compare_theme <- base_theme + theme(axis.title.y = element_blank())


# get thresholds for all predictors
thresholds <- getThresholdValues(pred_config, threshold_col = "alpha")
# filter merged data for CRISPR data to include and convert to format for bootstrapping
merged_training_bs <- merged_training %>% 
  filter(!is.na(direct_vs_indirect_negative)) %>%
  convertMergedForBootstrap(., pred_config = pred_config,
                            weight_col = "direct_vs_indirect_negative")

merged_heldout_bs <- merged_heldout %>% 
  filter(!is.na(direct_vs_indirect_negative)) %>%
  convertMergedForBootstrap(., pred_config = pred_config,
                            weight_col = "direct_vs_indirect_negative")


#### AUPRC ####
if (!plot_only) {
    ## Calculate delta performance between predictors for weighted AUPRC and precision
    weighted_delta_auprc_pred <- bootstrapDeltaPerformance(merged_heldout_bs, metric = "auprc", R = 1000, weighted = TRUE)
    weighted_delta_prec_pred <- bootstrapDeltaPerformance(merged_heldout_bs, metric = "precision", R = 1000, weighted = TRUE, thresholds = thresholds)
    write_tsv(weighted_delta_auprc_pred, file.path(out_dir, "heldout_weighted_AUPRC_delta_predictors.tsv"))
    write_tsv(weighted_delta_prec_pred, file.path(out_dir, "heldout_weighted_precision_delta_predictors.tsv"))

    # calculate bootstrapped un-weighted AUPRC for selected training and heldout data
    unweighted_auprc_training <- bootstrapPerformanceIntervals(merged_training_bs, metric = "auprc", R = 1000)
    unweighted_auprc_heldout <- bootstrapPerformanceIntervals(merged_heldout_bs, metric = "auprc", R = 1000)
    message("Calculated unweighted AUPRC CIs.")

    # calculate weighted AUPRC for selected training and heldout data
    weighted_auprc_training <- bootstrapPerformanceIntervals(merged_training_bs, metric = "auprc",
                                                            weighted = TRUE, R = 1000)
    weighted_auprc_heldout <- bootstrapPerformanceIntervals(merged_heldout_bs, metric = "auprc",
                                                            weighted = TRUE, R = 1000)
    message("Calculated weighted AUPRC CIs.")

    # combine into tables
    unweighted_auprc <- bind_rows(Training = unweighted_auprc_training, `Held-out` = unweighted_auprc_heldout,
        .id = "crispr_dataset")
    weighted_auprc <- bind_rows(Training = weighted_auprc_training, `Held-out` = weighted_auprc_heldout,
        .id = "crispr_dataset")

    # save tables to tsv
    print(file.path(out_dir, "unweighted_AUPRC.tsv"))
    print(file.path(out_dir, "weighted_AUPRC.tsv"))
    write_tsv(unweighted_auprc, file.path(out_dir, "unweighted_AUPRC.tsv"))
    write_tsv(weighted_auprc, file.path(out_dir, "weighted_AUPRC.tsv"))

    # compute bootstrapped un-weighted AUPRC differences between training and held-out data
    unweighted_delta_auprc <- bootstrapDeltaPerformanceDatasets(merged_training_bs, merged_heldout_bs,
                                                                metric = "auprc", R = 1000)

    # compute bootstrapped weighted AUPRC differences between training and held-out data
    weighted_delta_auprc <- bootstrapDeltaPerformanceDatasets(merged_training_bs, merged_heldout_bs,
                                                            metric = "auprc", weighted = TRUE, R = 1000)
    message("Calculated delta AUPRC significance.")                                 
    write_tsv(unweighted_delta_auprc, file.path(out_dir, "unweighted_delta_AUPRC.tsv"))
    write_tsv(weighted_delta_auprc, file.path(out_dir, "weighted_delta_AUPRC.tsv"))


} else {
    unweighted_auprc <- read_tsv(file.path(out_dir, "unweighted_AUPRC.tsv"), show_col_types = FALSE)
    weighted_auprc <- read_tsv(file.path(out_dir, "weighted_AUPRC.tsv"), show_col_types = FALSE)

    unweighted_delta_auprc <- read_tsv(file.path(out_dir, "unweighted_delta_AUPRC.tsv"), show_col_types = FALSE)
    weighted_delta_auprc <- read_tsv(file.path(out_dir, "weighted_delta_AUPRC.tsv"), show_col_types = FALSE)
}

# format files for plotting
unweighted_auprc <- unweighted_auprc %>% 
    mutate(pred_name_long = pred_key[id],
        fill_group = if_else(crispr_dataset == "Held-out", true = pred_name_long, false = "Training"))

# define pred order on unweighted auprc
pred_order <- pull(arrange(unweighted_auprc, full), pred_name_long) %>% unique()

unweighted_auprc <- mutate(unweighted_auprc,
    pred_name_long = factor(pred_name_long, levels = pred_order, ordered = TRUE))

unweighted_delta_auprc <- unweighted_auprc %>% 
    distinct(id, pred_name_long) %>% 
    inner_join(unweighted_delta_auprc, by = c("id")) %>% 
    mutate(pval_label = case_when(
                    pvalue < 0.001 ~ "***",
                    pvalue < 0.01 ~ "**",
                    pvalue < 0.05 ~ "*",
                    TRUE ~ "n.s."))

weighted_auprc <- weighted_auprc %>% 
    mutate(pred_name_long = pred_key[id],
        fill_group = if_else(crispr_dataset == "Held-out", true = pred_name_long, false = "Training"),
        pred_name_long = factor(pred_name_long, levels = pred_order, ordered = TRUE))

weighted_delta_auprc <- weighted_auprc %>% 
    distinct(id, pred_name_long) %>% 
    inner_join(weighted_delta_auprc, by = c("id")) %>% 
    mutate(pval_label = case_when(
                    pvalue < 0.001 ~ "***",
                    pvalue < 0.01 ~ "**",
                    pvalue < 0.05 ~ "*",
                    TRUE ~ "n.s."))

# plot auprc for each predictor both on the training and the held-out CRISPR data 
p_auprc <- ggplot(unweighted_auprc, aes(y = pred_name_long, x = full, color = pred_name_long,
                           fill = fill_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                position = position_dodge(0.9), width = 0.4, color = "black", alpha = 1) +
  geom_text(data = unweighted_delta_auprc,
        aes(y = pred_name_long, x = 0.85, label = pval_label),
        inherit.aes = FALSE, size = 3.2, hjust = 1, vjust = 0.5) +
  labs(title = "AUPRC", x = "AUPRC", subtitle = unweighted_title,
       alpha = "CRISPR\ndataset", fill = "Predictor") +
  scale_color_manual(values = pred_colors) +
  scale_fill_manual(values = fill_colors, na.value = NA) +
  scale_x_continuous(limits = c(0, 1)) +
  guides(color = "none") +
  theme_classic() + compare_theme

ggsave(filename = file.path(out_dir, "auprc_training_vs_heldout.pdf"), p_auprc,
       height = h, width = w)

p_auprc_weighted <- ggplot(weighted_auprc, aes(y = pred_name_long, x = full, color = pred_name_long,
                           fill = fill_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                position = position_dodge(0.9), width = 0.4, color = "black", alpha = 1) +
  geom_text(data = weighted_delta_auprc,
        aes(y = pred_name_long, x = 0.85, label = pval_label),
        inherit.aes = FALSE, size = 3.2, hjust = 1, vjust = 0.5) +
  labs(title = "Weighted AUPRC", x = "Weighted AUPRC", subtitle = weighted_title,
       alpha = "CRISPR\ndataset", fill = "Predictor") +
  scale_color_manual(values = pred_colors) +
  scale_fill_manual(values = fill_colors, na.value = NA) +
  scale_x_continuous(limits = c(0, 1)) +
  guides(color = "none") +
  theme_classic() + compare_theme

# save plot to pdf
ggsave(filename = file.path(out_dir, "weighted_auprc_training_vs_heldout.pdf"), p_auprc_weighted,
       height = h, width = w)
message("Plotted weighted AUPRC.")


#### PRECISION ####
if (!plot_only) {
    # calculate unweighted precision for selected training and heldout data
    unweighted_precision_training <- bootstrapPerformanceIntervals(merged_training_bs, metric = "precision", 
                                                                thresholds = thresholds, R = 1000)
    unweighted_precision_heldout <- bootstrapPerformanceIntervals(merged_heldout_bs, metric = "precision", 
                                                            thresholds = thresholds, R = 1000)
    message("Calculated unweighted precisions.")

    # calculate weighted precision for selected training and heldout data
    weighted_precision_training <- bootstrapPerformanceIntervals(merged_training_bs, metric = "precision", 
                                                            thresholds = thresholds, weighted = TRUE,
                                                            R = 1000)
    weighted_precision_heldout <- bootstrapPerformanceIntervals(merged_heldout_bs, metric = "precision", 
                                                            thresholds = thresholds, weighted = TRUE,
                                                            R = 1000)
    message("Calculated weighted precisions.")

    # combine into tables
    weighted_precision <- bind_rows(Training = weighted_precision_training, `Held-out` = weighted_precision_heldout,
                                .id = "crispr_dataset")
    unweighted_precision <- bind_rows(Training = unweighted_precision_training,
                                `Held-out` = unweighted_precision_heldout, .id = "crispr_dataset")
    write_tsv(unweighted_precision, file.path(out_dir, "unweighted_precision.tsv"))
    write_tsv(weighted_precision, file.path(out_dir, "weighted_precision.tsv"))

    # compute bootstrapped un-weighted precision differences between training and held-out data
    unweighted_delta_precision <- bootstrapDeltaPerformanceDatasets(merged_training_bs, merged_heldout_bs,
                                                                metric = "precision",
                                                                thresholds = thresholds, R = 1000)

    # compute bootstrapped weighted precision differences between training and held-out data
    weighted_delta_precision <- bootstrapDeltaPerformanceDatasets(merged_training_bs, merged_heldout_bs,
                                                            metric = "precision",
                                                            thresholds = thresholds, weighted = TRUE,
                                                            R = 1000)                                                        
    message("Calculated delta precision significances.")
    write_tsv(unweighted_delta_precision, file.path(out_dir, "unweighted_delta_precision.tsv"))
    write_tsv(weighted_delta_precision, file.path(out_dir, "weighted_delta_precision.tsv"))


} else {
    unweighted_precision <- read_tsv(file.path(out_dir, "unweighted_precision.tsv"), show_col_types = FALSE)
    weighted_precision <- read_tsv(file.path(out_dir, "weighted_precision.tsv"), show_col_types = FALSE)

    unweighted_delta_precision <- read_tsv(file.path(out_dir, "unweighted_delta_precision.tsv"), show_col_types = FALSE)
    weighted_delta_precision <- read_tsv(file.path(out_dir, "weighted_delta_precision.tsv"), show_col_types = FALSE)
}

# format inputs for plotting
weighted_precision <- weighted_precision %>%
    mutate(pred_name_long = pred_key[id],
        fill_group = if_else(crispr_dataset == "Held-out", true = pred_name_long, false = "Training"),
        pred_name_long = factor(pred_name_long, levels = pred_order, ordered = TRUE))

weighted_delta_precision <- weighted_precision %>% 
    distinct(id, pred_name_long) %>% 
    inner_join(weighted_delta_precision, by = c("id")) %>% 
    mutate(pval_label = case_when(
                    pvalue < 0.001 ~ "***",
                    pvalue < 0.01 ~ "**",
                    pvalue < 0.05 ~ "*",
                    TRUE ~ "n.s."))

unweighted_precision <- unweighted_precision %>%
    mutate(pred_name_long = pred_key[id],
        fill_group = if_else(crispr_dataset == "Held-out", true = pred_name_long, false = "Training"),
        pred_name_long = factor(pred_name_long, levels = pred_order, ordered = TRUE))

unweighted_delta_precision <- unweighted_precision %>% 
    distinct(id, pred_name_long) %>% 
    inner_join(unweighted_delta_precision, by = c("id")) %>% 
    mutate(pval_label = case_when(
                    pvalue < 0.001 ~ "***",
                    pvalue < 0.01 ~ "**",
                    pvalue < 0.05 ~ "*",
                    TRUE ~ "n.s."))
# plot
p_prec <- ggplot(unweighted_precision, aes(y = pred_name_long, x = full, color = pred_name_long,
                           fill = fill_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                position = position_dodge(0.9), width = 0.4, color = "black", alpha = 1) +
  geom_text(data = unweighted_delta_precision,
        aes(y = pred_name_long, x = 0.8, label = pval_label),
        inherit.aes = FALSE, size = 3.2, hjust = 1, vjust = 0.5) +
  labs(title = "Precision at threshold", x = "Precision at threshold", subtitle = unweighted_title,
       alpha = "CRISPR\ndataset", fill = "Predictor") +
  scale_color_manual(values = pred_colors) +
  scale_fill_manual(values = fill_colors, na.value = NA) +
  scale_x_continuous(limits = c(0, 1)) +
  guides(color = "none") +
  theme_classic() + compare_theme
ggsave(filename = file.path(out_dir, "precision_training_vs_heldout.pdf"), p_prec,
       height = h, width = w)
message("Plotted unweighted precisions.")

p_prec_weighted <- ggplot(weighted_precision, aes(y = pred_name_long, x = full, color = pred_name_long,
                           fill = fill_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                position = position_dodge(0.9), width = 0.4, color = "black", alpha = 1) +
  geom_text(data = weighted_delta_precision,
        aes(y = pred_name_long, x = 0.8, label = pval_label),
        inherit.aes = FALSE, size = 3.2, hjust = 1, vjust = 0.5) +
  labs(title = "Weighted precision at threshold", x = "Weighted precision at threshold", subtitle = weighted_title,
       alpha = "CRISPR\ndataset", fill = "Predictor") +
  scale_color_manual(values = pred_colors) +
  scale_fill_manual(values = fill_colors, na.value = NA) +
  scale_x_continuous(limits = c(0, 1)) +
  guides(color = "none") +
  theme_classic() + compare_theme
ggsave(filename = file.path(out_dir, "weighted_precision_training_vs_heldout.pdf"), p_prec_weighted,
       height = h, width = w)
message("Plotted weighted precisions.")

#### RECALL ####
if (!plot_only) {
    # calculate unweighted recall for selected training and heldout data
    unweighted_recall_training <- bootstrapPerformanceIntervals(merged_training_bs, metric = "recall", 
                                                                thresholds = thresholds, R = 1000)
    unweighted_recall_heldout <- bootstrapPerformanceIntervals(merged_heldout_bs, metric = "recall", 
                                                            thresholds = thresholds, R = 1000)
    message("Calculated unweighted recall CIs.")

    # calculate weighted recall for selected training and heldout data
    weighted_recall_training <- bootstrapPerformanceIntervals(merged_training_bs, metric = "recall", 
                                                            thresholds = thresholds, weighted = TRUE,
                                                            R = 1000)
    weighted_recall_heldout <- bootstrapPerformanceIntervals(merged_heldout_bs, metric = "recall", 
                                                            thresholds = thresholds, weighted = TRUE,
                                                            R = 1000)
    message("Calculated weighted recall CIs.")

    # combine into tables
    weighted_recall <- bind_rows(Training = weighted_recall_training, `Held-out` = weighted_recall_heldout,
                                .id = "crispr_dataset")
    unweighted_recall <- bind_rows(Training = unweighted_recall_training,
                                `Held-out` = unweighted_recall_heldout, .id = "crispr_dataset")
    write_tsv(unweighted_recall, file.path(out_dir, "unweighted_recall.tsv"))
    write_tsv(weighted_recall, file.path(out_dir, "weighted_recall.tsv"))

    # compute bootstrapped un-weighted AUPRC differences between training and held-out data
    unweighted_delta_recall <- bootstrapDeltaPerformanceDatasets(merged_training_bs, merged_heldout_bs,
                                                                metric = "recall",
                                                                thresholds = thresholds, R = 1000)

    # compute bootstrapped weighted AUPRC differences between training and held-out data
    weighted_delta_recall <- bootstrapDeltaPerformanceDatasets(merged_training_bs, merged_heldout_bs,
                                                            metric = "recall",
                                                            thresholds = thresholds, weighted = TRUE,
                                                            R = 1000)            
    message("Calculated delta recall significance.")
                                                                                                        
    write_tsv(unweighted_delta_recall, file.path(out_dir, "unweighted_delta_recall.tsv"))
    write_tsv(weighted_delta_recall, file.path(out_dir, "weighted_delta_recall.tsv"))

} else {
    weighted_recall <- read_tsv(file.path(out_dir, "weighted_recall.tsv"), show_col_types = FALSE)
    unweighted_recall <- read_tsv(file.path(out_dir, "unweighted_recall.tsv"), show_col_types = FALSE)

    unweighted_delta_recall <- read_tsv(file.path(out_dir, "unweighted_delta_recall.tsv"), show_col_types = FALSE)
    weighted_delta_recall <- read_tsv(file.path(out_dir, "weighted_delta_recall.tsv"), show_col_types = FALSE)
}

# format inputs for plotting
weighted_recall <- weighted_recall %>%
    mutate(pred_name_long = pred_key[id],
        fill_group = if_else(crispr_dataset == "Held-out", true = pred_name_long, false = "Training"),
        pred_name_long = factor(pred_name_long, levels = pred_order, ordered = TRUE))

weighted_delta_recall <- weighted_recall %>% 
    distinct(id, pred_name_long) %>% 
    inner_join(weighted_delta_recall, by = c("id")) %>% 
    mutate(pval_label = case_when(
                    pvalue < 0.001 ~ "***",
                    pvalue < 0.01 ~ "**",
                    pvalue < 0.05 ~ "*",
                    TRUE ~ "n.s."))

unweighted_recall <- unweighted_recall %>% mutate(pred_name_long = pred_key[id],
    fill_group = if_else(crispr_dataset == "Held-out", true = pred_name_long, false = "Training"),
    pred_name_long = factor(pred_name_long, levels = pred_order, ordered = TRUE))

unweighted_delta_recall <- unweighted_recall %>% 
    distinct(id, pred_name_long) %>% 
    inner_join(unweighted_delta_recall, by = c("id")) %>% 
    mutate(pval_label = case_when(
                    pvalue < 0.001 ~ "***",
                    pvalue < 0.01 ~ "**",
                    pvalue < 0.05 ~ "*",
                    TRUE ~ "n.s."))


# plot
p_recall <- ggplot(unweighted_recall, aes(y = pred_name_long, x = full, color = pred_name_long,
                           fill = fill_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                position = position_dodge(0.9), width = 0.4, color = "black", alpha = 1) +
  geom_text(data = unweighted_delta_recall,
        aes(y = pred_name_long, x = 0.8, label = pval_label),
        inherit.aes = FALSE, size = 3.2, hjust = 1, vjust = 0.5) +
  labs(title = "Recall at threshold", x = "Recall at threshold", subtitle = unweighted_title,
       alpha = "CRISPR\ndataset", fill = "Predictor") +
  scale_color_manual(values = pred_colors) +
  scale_fill_manual(values = fill_colors, na.value = NA) +
  scale_x_continuous(limits = c(0, 1)) +
  guides(color = "none") +
  theme_classic() + compare_theme
ggsave(filename = file.path(out_dir, "recall_training_vs_heldout.pdf"), p_recall,
       height = h, width = w)
message("Plotted unweighted recalls.")

p_recall_weighted <- ggplot(weighted_recall, aes(y = pred_name_long, x = full, color = pred_name_long,
                           fill = fill_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                position = position_dodge(0.9), width = 0.4, color = "black", alpha = 1) +
  geom_text(data = weighted_delta_recall,
        aes(y = pred_name_long, x = 0.8, label = pval_label),
        inherit.aes = FALSE, size = 3.2, hjust = 1, vjust = 0.5) +
  labs(title = "Weighted recall at threshold", x = "Weighted recall at threshold", subtitle = weighted_title,
       alpha = "CRISPR\ndataset", fill = "Predictor") +
  scale_color_manual(values = pred_colors) +
  scale_fill_manual(values = fill_colors, na.value = NA) +
  scale_x_continuous(limits = c(0, 1)) +
  guides(color = "none") +
  theme_classic() + compare_theme
ggsave(filename = file.path(out_dir, "weighted_recall_training_vs_heldout.pdf"), p_recall_weighted,
       height = h, width = w)
message("Plotted weighted recalls.")

### final plot of auprc and precision just on held-out
weighted_auprc <- read_tsv(file.path(out_dir, "weighted_AUPRC.tsv"), show_col_types = FALSE) %>%
    filter(crispr_dataset != "Training") %>% 
    mutate(pred_name_long = pred_key[id],
                pred_name_long = factor(pred_name_long, levels = rev(pred_order), ordered = TRUE))
weighted_precision <- read_tsv(file.path(out_dir, "weighted_precision.tsv"), show_col_types = FALSE) %>%
    filter(crispr_dataset != "Training") %>% 
    mutate(pred_name_long = pred_key[id],
        pred_name_long = factor(pred_name_long, levels = rev(pred_order), ordered = TRUE))

p2_auprc <- ggplot(weighted_auprc, aes(x = pred_name_long, y = full, fill = pred_name_long)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(0.9), width = 0.4, color = "black", alpha = 1) +
  labs(title = "Weighted AUPRC", y = "Weighted AUPRC", fill = "Predictor", subtitle = weighted_heldout_title) +
  scale_fill_manual(values = pred_colors, na.value = NA) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic() + base_theme +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
    axis.text = element_text(color = "#000000", size = 13),
    legend.position = "none")

p2_prec <- ggplot(weighted_precision, aes(x = pred_name_long, y = full, fill = pred_name_long)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(0.9), width = 0.4, color = "black", alpha = 1) +
  labs(title = "Weighted precision", y = "Weighted precision at threshold", fill = "Predictor", subtitle = weighted_heldout_title) +
  scale_fill_manual(values = pred_colors, na.value = NA) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic() + base_theme +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text = element_text(color = "#000000", size = 13),
    legend.position = "none")

grid <- cowplot::plot_grid(p2_auprc, p2_prec, nrow = 1, rel_widths = c(1, 1), align = "hv")
ggsave2(file.path(out_dir, "weighted_auprc_precision_heldout_only.pdf"), height = 5, width = 7)
