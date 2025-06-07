library(tidyverse)
library(cowplot)

proj_dir <- "/oak/stanford/groups/engreitz/Users/sheth/CRISPR_comparison_v3/CRISPR_comparison"
source(file.path(proj_dir, "workflow/scripts/crisprComparisonLoadInputData.R"))
source(file.path(proj_dir, "workflow/scripts/crisprComparisonPlotFunctions.R"))
source(file.path(proj_dir, "workflow/scripts/bootstrap/crisprComparisonBootstrapFunctions_weighted_trapAUC.R"))

## file paths
merged_heldout_crispr_file <- "/oak/stanford/groups/engreitz/Projects/scE2G/data/CRISPR_benchmarking/held_out_250602/expt_pred_merged_annot.txt.gz"
merged_training_crispr_file <- "/oak/stanford/groups/engreitz/Projects/scE2G/data/CRISPR_benchmarking/sc.250425_2.K562.Xu.cv/expt_pred_merged_annot.txt.gz"
pred_config_file <- file.path(proj_dir, "resources", "pred_config", "pred_config_scE2G_heldout.tsv")
correct_heldout_annotations_file <- "/oak/stanford/groups/engreitz/Users/agschwin/distal_regulation_paper/CRISPR_indirect_effects/results/annotated_crispr_data/EPCrisprBenchmark_combined_heldout_element_classes_direct_effects.tsv.gz"
correct_training_annotations_file <- "/oak/stanford/groups/engreitz/Users/agschwin/distal_regulation_paper/CRISPR_indirect_effects/results/annotated_crispr_data/EPCrisprBenchmark_combined_training_element_classes_direct_effects.tsv.gz"


## load input files
merged_heldout <- read_tsv(merged_heldout_crispr_file, show_col_types = FALSE) %>%
    select(-any_of(c("element_category_with_dnase", "element_category_simple",	"element_category", "distance", "distanceToTSS", "Dataset")))

merged_training <- read_tsv(merged_training_crispr_file, show_col_types = FALSE) %>% 
    select(-any_of(c("element_category_with_dnase", "element_category_simple",	"element_category", "distance", "distanceToTSS", "Dataset"))) %>%
    mutate(pred_uid = ifelse(pred_uid == "scE2G_ATAC.E2G.Score.cv.qnorm", "scE2G_ATAC.E2G.Score.qnorm", pred_uid),
        pred_uid = ifelse(pred_uid == "scE2G_multiome.E2G.Score.cv.qnorm", "scE2G_multiome.E2G.Score.qnorm", pred_uid),
        pred_col = ifelse(pred_col == "E2G.Score.cv.qnorm", "E2G.Score.qnorm", pred_col))

pred_config <- importPredConfig(pred_config_file, expr = FALSE)
pred_key <- setNames(pred_config$pred_name_long, pred_config$pred_uid)

## add new chromatin categories and direct effect rate
## filter to specific chromatin categoreis and/or effect sizes, distance to tss
out_dir <- file.path(proj_dir, "workflow", "results", "v2_scE2G_heldout_analysis"); dir.create(out_dir, showWarnings = FALSE)
out_dir <- file.path(out_dir, "set_CTCF_H3K27me3_false_filter_low2K27ac_ES5pct");  dir.create(out_dir, showWarnings = FALSE)

correct_annot_heldout <- read_tsv(correct_heldout_annotations_file) %>%
    select(chrom, chromStart, chromEnd, measuredGeneSymbol, ExperimentCellType = CellType, element_category, direct_vs_indirect_negative, distanceToTSS, Dataset)
correct_annot_training <- read_tsv(correct_training_annotations_file) %>%
    select(chrom, chromStart, chromEnd, measuredGeneSymbol, ExperimentCellType = CellType, element_category, direct_vs_indirect_negative, distanceToTSS, Dataset)

merged_heldout <- left_join(merged_heldout, correct_annot_heldout, by = c("chrom", "chromStart", "chromEnd", "measuredGeneSymbol", "ExperimentCellType")) %>% 
    filter(distanceToTSS < 1000000) %>% 
    mutate(Regulated = ifelse(element_category %in% c("CTCF overlap", "H3K27me3 overlap"), FALSE, Regulated)) %>% 
    filter(!(Regulated & abs(EffectSize) <= 0.05)) %>% 
    filter(!(Regulated & element_category == "H3K27ac low"))
    #filter(ExperimentCellType != "K562")

merged_training <- left_join(merged_training, correct_annot_training, by = c("chrom", "chromStart", "chromEnd", "measuredGeneSymbol", "ExperimentCellType")) %>% 
    filter(distanceToTSS < 1000000) %>% 
    mutate(Regulated = ifelse(element_category %in% c("CTCF overlap", "H3K27me3 overlap"), FALSE, Regulated)) %>% 
    filter(!(Regulated & abs(EffectSize) <= 0.05)) %>% 
    filter(!(Regulated & element_category == "H3K27ac low"))

heldout_summary <- merged_heldout %>%
    select(chrom, chromStart, chromEnd, measuredGeneSymbol, Regulated, ExperimentCellType, Dataset) %>% 
    distinct()
training_summary <- merged_training %>% 
    select(chrom, chromStart, chromEnd, measuredGeneSymbol, Regulated, ExperimentCellType, Dataset) %>% 
    distinct()
combined_summary <- rbind(heldout_summary, training_summary) %>% 
    group_by(ExperimentCellType, Dataset) %>% 
    summarize(n_total = n(), n_regulated = sum(Regulated))
write_tsv(combined_summary, file.path(out_dir, "crispr_dataset_summary.tsv"))

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

### Calculate delta performance between predictors for weighted AUPRC and precision
weighted_delta_auprc_pred <- bootstrapDeltaPerformance(merged_heldout_bs, metric = "auprc", R = 1000, weighted = TRUE)
weighted_delta_prec_pred <- bootstrapDeltaPerformance(merged_heldout_bs, metric = "precision", R = 1000, weighted = TRUE, thresholds = thresholds)
write_tsv(weighted_delta_auprc_pred, file.path(out_dir, "heldout_weighted_AUPRC_delta_predictors.tsv"))
write_tsv(weighted_delta_prec_pred, file.path(out_dir, "heldout_weighted_precision_delta_predictors.tsv"))

### Weighted performance 
## AUPRC
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
unweighted_auprc <- bind_rows(Training = unweighted_auprc_training,
                              `Held-out` = unweighted_auprc_heldout, .id = "crispr_dataset")
weighted_auprc <- bind_rows(Training = weighted_auprc_training, `Held-out` = weighted_auprc_heldout,
                            .id = "crispr_dataset") %>% 

# save tables to pdf
print(file.path(out_dir, "unweighted_AUPRC.tsv"))
print(file.path(out_dir, "weighted_AUPRC.tsv"))
write_tsv(unweighted_auprc, file.path(out_dir, "unweighted_AUPRC.tsv"))
write_tsv(weighted_auprc, file.path(out_dir, "weighted_AUPRC.tsv"))

unweighted_auprc <- read_tsv(file.path(out_dir, "unweighted_AUPRC.tsv"))
weighted_auprc <- read_tsv(file.path(out_dir, "weighted_AUPRC.tsv"))

unweighted_auprc <- unweighted_auprc %>% 
    mutate(pred_name_long = pred_key[id],
        fill_group = if_else(crispr_dataset == "Held-out", true = pred_name_long, false = "Training"))

# define pred order on unweighted auprc
pred_order <- pull(arrange(unweighted_auprc, full), pred_name_long) %>% unique()

unweighted_auprc <- mutate(unweighted_auprc,
    pred_name_long = factor(pred_name_long, levels = pred_order, ordered = TRUE))

weighted_auprc <- weighted_auprc %>% 
    mutate(pred_name_long = pred_key[id],
        fill_group = if_else(crispr_dataset == "Held-out", true = pred_name_long, false = "Training"),
        pred_name_long = factor(pred_name_long, levels = pred_order, ordered = TRUE))

# plot recall for each predictor both on the training and the held-out CRISPR data 
p_auprc <- ggplot(unweighted_auprc, aes(y = pred_name_long, x = full, color = pred_name_long,
                           fill = fill_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                position = position_dodge(0.9), width = 0.4, color = "black", alpha = 1) +
  labs(title = "AUPRC", x = "AUPRC",
       alpha = "CRISPR\ndataset", fill = "Predictor") +
  scale_color_manual(values = pred_colors) +
  scale_fill_manual(values = fill_colors, na.value = NA) +
  scale_x_continuous(limits = c(0, 1)) +
  guides(color = "none") +
  theme_classic() +
  theme(axis.title.y = element_blank(), text = element_text(size = 13))
ggsave(filename = file.path(out_dir, "auprc_training_vs_heldout.pdf"), p_auprc,
       height = 8, width = 12)

p_auprc_weighted <- ggplot(weighted_auprc, aes(y = pred_name_long, x = full, color = pred_name_long,
                           fill = fill_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                position = position_dodge(0.9), width = 0.4, color = "black", alpha = 1) +
  labs(title = "Weighted AUPRC", x = "Weighted AUPRC",
       alpha = "CRISPR\ndataset", fill = "Predictor") +
  scale_color_manual(values = pred_colors) +
  scale_fill_manual(values = fill_colors, na.value = NA) +
  scale_x_continuous(limits = c(0, 1)) +
  guides(color = "none") +
  theme_classic() +
  theme(axis.title.y = element_blank(), text = element_text(size = 13))

# save plot to pdf
ggsave(filename = file.path(out_dir, "weighted_auprc_training_vs_heldout.pdf"), p_auprc_weighted,
       height = 8, width = 12)
message("Plotted weighted AUPRC.")

# compute bootstrapped un-weighted AUPRC differences between training and held-out data
unweighted_delta_auprc <- bootstrapDeltaPerformanceDatasets(merged_training_bs, merged_heldout_bs,
                                                            metric = "auprc", R = 1000)

# compute bootstrapped weighted AUPRC differences between training and held-out data
weighted_delta_auprc <- bootstrapDeltaPerformanceDatasets(merged_training_bs, merged_heldout_bs,
                                                          metric = "auprc", weighted = TRUE, R = 1000)
message("Calculated delta AUPRC significance.")                                 
write_tsv(unweighted_delta_auprc, file.path(out_dir, "unweighted_delta_AUPRC.tsv"))
write_tsv(weighted_delta_auprc, file.path(out_dir, "weighted_delta_AUPRC.tsv"))

## Recall
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

# plot recalls
weighted_recall <- weighted_recall %>%
    mutate(pred_name_long = pred_key[id],
        fill_group = if_else(crispr_dataset == "Held-out", true = pred_name_long, false = "Training"),
        pred_name_long = factor(pred_name_long, levels = pred_order, ordered = TRUE))
unweighted_recall <- unweighted_recall %>% mutate(pred_name_long = pred_key[id],
    fill_group = if_else(crispr_dataset == "Held-out", true = pred_name_long, false = "Training"),
    pred_name_long = factor(pred_name_long, levels = pred_order, ordered = TRUE))

p_recall <- ggplot(unweighted_recall, aes(y = pred_name_long, x = full, color = pred_name_long,
                           fill = fill_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                position = position_dodge(0.9), width = 0.4, color = "black", alpha = 1) +
  labs(title = "Recall at threshold", x = "Recall at threshold",
       alpha = "CRISPR\ndataset", fill = "Predictor") +
  scale_color_manual(values = pred_colors) +
  scale_fill_manual(values = fill_colors, na.value = NA) +
  scale_x_continuous(limits = c(0, 1)) +
  guides(color = "none") +
  theme_classic() +
  theme(axis.title.y = element_blank(), text = element_text(size = 13))
ggsave(filename = file.path(out_dir, "recall_training_vs_heldout.pdf"), p_recall,
       height = 8, width = 12)
message("Plotted unweighted recalls.")

p_recall_weighted <- ggplot(weighted_recall, aes(y = pred_name_long, x = full, color = pred_name_long,
                           fill = fill_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                position = position_dodge(0.9), width = 0.4, color = "black", alpha = 1) +
  labs(title = "Weighted recall at threshold", x = "Weighted recall at threshold",
       alpha = "CRISPR\ndataset", fill = "Predictor") +
  scale_color_manual(values = pred_colors) +
  scale_fill_manual(values = fill_colors, na.value = NA) +
  scale_x_continuous(limits = c(0, 1)) +
  guides(color = "none") +
  theme_classic() +
  theme(axis.title.y = element_blank(), text = element_text(size = 13))
ggsave(filename = file.path(out_dir, "weighted_recall_training_vs_heldout.pdf"), p_recall_weighted,
       height = 8, width = 12)
message("Plotted weighted recalls.")

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

## Precision
# calculate unweighted recall for selected training and heldout data
unweighted_precision_training <- bootstrapPerformanceIntervals(merged_training_bs, metric = "precision", 
                                                            thresholds = thresholds, R = 1000)
unweighted_precision_heldout <- bootstrapPerformanceIntervals(merged_heldout_bs, metric = "precision", 
                                                           thresholds = thresholds, R = 1000)
message("Calculated unweighted precisions.")

# calculate weighted recall for selected training and heldout data
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

unweighted_precision <- read_tsv(file.path(out_dir, "unweighted_precision.tsv"))
weighted_precision <- read_tsv(file.path(out_dir, "weighted_precision.tsv"))

# plot recalls
weighted_precision <- weighted_precision %>%
    mutate(pred_name_long = pred_key[id],
        fill_group = if_else(crispr_dataset == "Held-out", true = pred_name_long, false = "Training"),
        pred_name_long = factor(pred_name_long, levels = pred_order, ordered = TRUE))
unweighted_precision <- unweighted_precision %>%
    mutate(pred_name_long = pred_key[id],
        fill_group = if_else(crispr_dataset == "Held-out", true = pred_name_long, false = "Training"),
        pred_name_long = factor(pred_name_long, levels = pred_order, ordered = TRUE))

p_prec <- ggplot(unweighted_precision, aes(y = pred_name_long, x = full, color = pred_name_long,
                           fill = fill_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                position = position_dodge(0.9), width = 0.4, color = "black", alpha = 1) +
  labs(title = "Precision at threshold", x = "Precision at threshold",
       alpha = "CRISPR\ndataset", fill = "Predictor") +
  scale_color_manual(values = pred_colors) +
  scale_fill_manual(values = fill_colors, na.value = NA) +
  scale_x_continuous(limits = c(0, 1)) +
  guides(color = "none") +
  theme_classic() +
  theme(axis.title.y = element_blank(), text = element_text(size = 13))
ggsave(filename = file.path(out_dir, "precision_training_vs_heldout.pdf"), p_prec,
       height = 8, width = 12)
message("Plotted unweighted precisions.")

p_prec_weighted <- ggplot(weighted_precision, aes(y = pred_name_long, x = full, color = pred_name_long,
                           fill = fill_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                position = position_dodge(0.9), width = 0.4, color = "black", alpha = 1) +
  labs(title = "Weighted precision at threshold", x = "Weighted precision at threshold",
       alpha = "CRISPR\ndataset", fill = "Predictor") +
  scale_color_manual(values = pred_colors) +
  scale_fill_manual(values = fill_colors, na.value = NA) +
  scale_x_continuous(limits = c(0, 1)) +
  guides(color = "none") +
  theme_classic() +
  theme(axis.title.y = element_blank(), text = element_text(size = 13))
ggsave(filename = file.path(out_dir, "weighted_precision_training_vs_heldout.pdf"), p_prec_weighted,
       height = 8, width = 12)
message("Plotted weighted precisions.")

# compute bootstrapped un-weighted AUPRC differences between training and held-out data
unweighted_delta_precision <- bootstrapDeltaPerformanceDatasets(merged_training_bs, merged_heldout_bs,
                                                             metric = "precision",
                                                             thresholds = thresholds, R = 1000)

# compute bootstrapped weighted AUPRC differences between training and held-out data
weighted_delta_precision <- bootstrapDeltaPerformanceDatasets(merged_training_bs, merged_heldout_bs,
                                                           metric = "precision",
                                                           thresholds = thresholds, weighted = TRUE,
                                                           R = 1000)                                                        
message("Calculated delta precision significances.")
write_tsv(unweighted_delta_precision, file.path(out_dir, "unweighted_delta_precision.tsv"))
write_tsv(weighted_delta_precision, file.path(out_dir, "weighted_delta_precision.tsv"))

### final plot of auprc and precision just on held-out
weighted_auprc <- read_tsv(file.path(out_dir, "weighted_AUPRC.tsv")) %>%
    filter(crispr_dataset != "Training") %>% 
    mutate(pred_name_long = pred_key[id],
                pred_name_long = factor(pred_name_long, levels = rev(pred_order), ordered = TRUE))
weighted_precision <- read_tsv(file.path(out_dir, "weighted_precision.tsv")) %>%
    filter(crispr_dataset != "Training") %>% 
    mutate(pred_name_long = pred_key[id],
        pred_name_long = factor(pred_name_long, levels = rev(pred_order), ordered = TRUE))

p2_auprc <- ggplot(weighted_auprc, aes(x = pred_name_long, y = full, fill = pred_name_long)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(0.9), width = 0.4, color = "black", alpha = 1) +
  labs(title = "Weighted AUPRC", y = "Weighted AUPRC", fill = "Predictor") +
  scale_fill_manual(values = pred_colors, na.value = NA) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text = element_text(color = "#000000", size = 13),
    legend.position = "none")

p2_prec <- ggplot(weighted_precision, aes(x = pred_name_long, y = full, fill = pred_name_long)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(0.9), width = 0.4, color = "black", alpha = 1) +
  labs(title = "Weighted precision", y = "Weighted precision at threshold", fill = "Predictor") +
  scale_fill_manual(values = pred_colors, na.value = NA) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text = element_text(color = "#000000", size = 13),
    legend.position = "none")

grid <- cowplot::plot_grid(p2_auprc, p2_prec, nrow = 1, rel_widths = c(1, 1), align = "hv")
ggsave2(file.path(out_dir, "weighted_auprc_precision_heldout_only.pdf"), height = 5, width = 7)
