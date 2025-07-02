suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(tidyverse)
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(egg)
})

map_values_to_colors <- function(values, gradient, color_lims, na_color) {
  
  # create a color ramp function based on the color gradient
  gradient <- colorRamp(gradient)
  
  # scale the values to be between 0 and 1 based on color limits, using scales::rescale
  scaled_values <- scales::rescale(values, to = c(0, 1), from = color_lims)
  
  # squish out-of-bound values to the nearest limit
  scaled_values <- scales::squish(scaled_values, range = c(0, 1))
  
  # generate colors based on scaled values
  colors <- rgb(gradient(scaled_values), maxColorValue = 255)
  
  # handle NA values by assigning them the specified NA color
  colors[is.na(values)] <- na_color
  
  return(colors)
}

plot_single_correlation <- function(correlation_col, res_id, datasets_include, clusters_exclude) {
	input_file <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0214_new_global_properties/enhancer_similarity/correlation_across_clusters.tsv"
	sample_key_file <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0214_new_global_properties/config/pred_sample_key.tsv"
	# biosample	dataset	biosample_name	ext_biosample_name	pred_file	ATAC_bw
	
	out_dir_base <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0214_new_global_properties/enhancer_similarity"
	out_dir <- file.path(out_dir_base, res_id); dir.create(out_dir, showWarnings = FALSE)
	correlation_col <- correlation_cols[i]

	# format input data
	sample_key <- fread(sample_key_file) %>%
		dplyr::select(biosample, dataset, biosample_name)
	sample_key_A <- sample_key %>% setNames(c("biosampleA", "datasetA", "biosampleNameA"))
	sample_key_B <- sample_key %>% setNames(c("biosampleB", "datasetB", "biosampleNameB"))

	df <- fread(input_file) %>%
		distinct() %>%
		dplyr::filter(!(biosampleA %in% clusters_exclude), !(biosampleB %in% clusters_exclude)) %>%
		mutate(jaccard = nSharedPredAwB / (nTotalPredA + nTotalPredB - nSharedPredAwB)) %>%
		dplyr::select(biosampleA, biosampleB, !!sym(correlation_col))
	df_copy <- setNames(df, c("biosampleB", "biosampleA", correlation_col))
	df <- rbind(df, df_copy) %>%
		left_join(sample_key_A, by="biosampleA") %>%
		left_join(sample_key_B, by="biosampleB") %>%
		distinct() %>% 
		filter(datasetA %in% datasets_include, datasetB %in% datasets_include)

	df[[correlation_col]] <- as.numeric(df[[correlation_col]])
	print(head(df))
	#fwrite(df, out_table, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

	# cluster biosamples by similarity; save order
	matrix <- df %>% dplyr::select(biosampleNameA, biosampleNameB, !!sym(correlation_col)) %>%
		pivot_wider(names_from = biosampleNameA, values_from = !!sym(correlation_col)) %>%
		column_to_rownames("biosampleNameB")
	dist_mtx <- dist(matrix)
	print(dist_mtx)
	clust <- dist_mtx %>% replace_na(0) %>% hclust(method = "average")
	ordered_biosample_names <- rownames(matrix)[clust$order]
	saveRDS(ordered_biosample_names, file.path(out_dir, paste0("ordered_biosample_names_by_", correlation_col, ".rds")))

	df <- mutate(df,
		biosampleNameA =  factor(biosampleNameA, levels = ordered_biosample_names, ordered = TRUE),
		biosampleNameB = factor(biosampleNameB, levels = ordered_biosample_names, ordered = TRUE))

	# plotting params
	#gradient <- c("#f6eff7","#bdc9e1", "#67a9cf","#1c9099", "#016c59") # white - teal
	teals <- c("#cae5ee", "#96ced3",	"#49bcbc", "#0096a0", "#006479", "#003648") # monochrome teal
	blues <- c("#c5e5fb", "#9bcae9", "#5496ce", "#006eae", "#00488d", "#002359")
	blue_red <- c("#002359", "#006eae", "#9bcae9",  "#e5e5e9",  "#e9a0a5", "#c5373d", "#730c0d") # blue to red
	reds <- c("#f6ceca", "#e9a0a5", "#dc6464", "#c5373d", "#9b241c", "#730c0d") # red
	purples <- c("#e9d3ea", "#d3a9ce", "#b778b3", "#a64791", "#792374", "#430b4e")

	if (correlation_col %in% c("Pearson", "jaccard")) {
		if (correlation_col == "jaccard") {gradient <- reds} else { gradient <- blues}
		pct_lims <- c(0, 0.6)
	} else {
		gradient <- blue_red
		pct_lims <- c(-1, 1)
	}

	na_color <- "#ffffff" # real white

	ds_cols <- c("#e96a00", "#e9c54e", "#5496ce", "#00488d", "#5eb342", "#1c6e2b", "#9b241c", "#002359", "#0e3716", "#730c0d")
	names(ds_cols) = c("K562", "Other", "BMMC5", "BMMC22", "PBMC5", "PBMC9", "Islets", "BMMC_all", "PBMC_all", "Islets_all")
	ds_cols <- ds_cols[names(ds_cols) %in% datasets_include]

	# heatmap
	hm <- ggplot(df, aes(x=biosampleNameB, y=biosampleNameA, fill = !!sym(correlation_col))) +
		geom_tile() +
		scale_fill_gradientn(colors=gradient, oob=scales::squish, na.value=na_color, limits=pct_lims, name = correlation_col) + #name="% shared\nenhancers") +
		xlab("") + ylab("") +
		coord_fixed(1) + 
		theme_minimal() + theme(axis.text = element_text(size = 6), axis.title = element_blank(), axis.text.x = element_text(angle=60, hjust=1),
			legend.position='top',  legend.direction='horizontal', legend.text=element_text(size=7), legend.title=element_text(size=7))

	# dataset bar
	sample_key <- dplyr::filter(sample_key, dataset %in% datasets_include, !(biosample %in% clusters_exclude)) %>%
		mutate(biosample_name = factor(biosample_name, levels = ordered_biosample_names, ordered = TRUE))
	print(sample_key)
	print(sample_key$biosample_name)
		#dataset = factor(dataset, levels = c("BMMC5", "BMMC22", "BMMC_all", "PBMC5", "PBMC9", "PBMC_all", "Islets", "Islets_all", "K562", "GM12878"), ordered = TRUE)) %>%
		#dplyr::filter(dataset %in% datasets_include)
	ds <- ggplot(sample_key, aes(x = 1, y = biosample_name, fill = dataset)) +
		geom_tile() +
		scale_fill_manual(values = ds_cols, name = "Dataset") +
		theme_minimal() + theme(axis.text = element_blank(), axis.title = element_blank(), axis.text.x = element_blank(),
			legend.position='right',  legend.direction='vertical', legend.text=element_text(size=7), legend.title=element_text(size=7))

		# combine and save
		assembled <- egg::ggarrange(hm, ds, nrow=1, ncol=2, widths=c(8, 0.3))
		ggsave(file.path(out_dir, paste0(correlation_col, "_heatmap.pdf")), assembled, width = 7, height = 6)
}

plot_double_correlation <- function(res_id, datasets_include, clusters_exclude) {
	input_file <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0214_new_global_properties/enhancer_similarity/correlation_across_clusters.tsv"
	sample_key_file <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0214_new_global_properties/config/pred_sample_key.tsv"	
	out_dir_base <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0214_new_global_properties/enhancer_similarity"
	out_dir <- file.path(out_dir_base, res_id); dir.create(out_dir, showWarnings = FALSE)

	# format input data
	sample_key <- fread(sample_key_file) %>%
		dplyr::select(biosample, dataset, biosample_name)
	sample_key_A <- sample_key %>% setNames(c("biosampleA", "datasetA", "biosampleNameA"))
	sample_key_B <- sample_key %>% setNames(c("biosampleB", "datasetB", "biosampleNameB"))

	df <- fread(input_file) %>%
		distinct() %>%
		dplyr::filter(!(biosampleA %in% clusters_exclude), !(biosampleB %in% clusters_exclude)) %>%
		mutate(jaccard = nSharedPredAwB / (nTotalPredA + nTotalPredB - nSharedPredAwB)) %>%
		dplyr::select(biosampleA, biosampleB, jaccard, Pearson)
	df_copy <- setNames(df, c("biosampleB", "biosampleA", "jaccard", "Pearson"))
	df <- rbind(df, df_copy) %>%
		left_join(sample_key_A, by="biosampleA") %>%
		left_join(sample_key_B, by="biosampleB") %>%
		mutate(product = jaccard * Pearson) %>%
		distinct() %>% 
		filter(datasetA %in% datasets_include, datasetB %in% datasets_include)

	

	# cluster biosamples by similarity 
	matrix <- df %>% dplyr::select(biosampleNameA, biosampleNameB, product) %>%
		pivot_wider(names_from = biosampleNameA, values_from = product) %>%
		column_to_rownames("biosampleNameB")
	dist <- dist(matrix)
	clust <- dist %>% replace_na(0) %>% hclust(method = "average")
	ordered_biosample_names <- rownames(matrix)[clust$order]
	saveRDS(ordered_biosample_names, file.path(out_dir, paste0("ordered_biosample_names_Pearson_jaccard.rds")))

	df_unique <- df %>% #dplyr::filter(df, biosampleA != biosampleB) %>%
		mutate(biosampleNameA =  factor(biosampleNameA, levels = ordered_biosample_names, ordered = TRUE),
			biosampleNameB = factor(biosampleNameB, levels = ordered_biosample_names, ordered = TRUE),
			is_upper = as.numeric(biosampleNameA) > as.numeric(biosampleNameB),
			is_diag = (biosampleNameA == biosampleNameB))

	# plotting params
	#gradient <- c("#f6eff7","#bdc9e1", "#67a9cf","#1c9099", "#016c59") # white - teal
	teals <- c("#cae5ee", "#96ced3",	"#49bcbc", "#0096a0", "#006479", "#003648") # monochrome teal
	blues <- c("#c5e5fb", "#9bcae9", "#5496ce", "#006eae", "#00488d", "#002359")
	purples <- c("#e9d3ea", "#d3a9ce", "#b778b3", "#a64791", "#792374", "#430b4e")
	blue_red <- c("#002359", "#006eae", "#9bcae9",  "#e5e5e9",  "#e9a0a5", "#c5373d", "#730c0d") # blue to red
	reds <- c("#f6ceca", "#e9a0a5", "#dc6464", "#c5373d", "#9b241c", "#730c0d") # red

	color_lims <- c(0, 0.6)
	jaccard_cols <- blues
	cor_cols <- reds
	na_color <- "#ffffff" # real white
	greys <- c("#96a0b3", "#6e788d", "#1c2a43")

	df_unique$fill_color <- map_values_to_colors(df_unique$Pearson, cor_cols, color_lims, na_color)
	df_unique$fill_color[df_unique$is_upper] <- map_values_to_colors(df_unique$jaccard[df_unique$is_upper],
		jaccard_cols, color_lims, na_color)
	df_unique$fill_color[df_unique$is_diag] <- greys[3]

	 ds_cols <- c("#e96a00", "#e9c54e", "#5496ce", "#00488d", "#5eb342", "#1c6e2b", "#9b241c", "#002359", "#0e3716", "#730c0d")
	 names(ds_cols) = c("K562", "Other", "BMMC5", "BMMC22", "PBMC5", "PBMC9", "Islets", "BMMC_all", "PBMC_all", "Islets_all")
	ds_cols <- ds_cols[names(ds_cols) %in% datasets_include]

	# heatmap
	hm <- ggplot(df_unique, aes(x = biosampleNameA, y = biosampleNameB, fill = fill_color)) +
		geom_tile() +
		scale_fill_identity() +
		xlab("") + ylab("") +
		scale_x_discrete(limits = ordered_biosample_names) +
		scale_y_discrete(limits = ordered_biosample_names) +
		coord_fixed(1) + 
		theme_minimal() + theme(axis.text = element_text(size = 6.8), axis.title = element_blank(), axis.text.x = element_text(angle=60, hjust=1),
			legend.position='none',  legend.direction='horizontal', legend.text=element_text(size=7), legend.title=element_text(size=7))

	# dataset bar
	sample_key <- dplyr::filter(sample_key, dataset %in% datasets_include, !(biosample %in% clusters_exclude)) %>%
		mutate(biosample_name = factor(biosample_name, levels = ordered_biosample_names, ordered = TRUE))
		#dataset = factor(dataset, levels = c("BMMC5", "BMMC22", "BMMC_all", "PBMC5", "PBMC9", "PBMC_all", "Islets", "Islets_all", "K562", "GM12878"), ordered = TRUE)) %>%
		#dplyr::filter(dataset %in% datasets_include)
	ds <- ggplot(sample_key, aes(x = 1, y = biosample_name, fill = dataset)) +
		geom_tile() +
		scale_fill_manual(values = ds_cols, name = "Dataset") +
		theme_minimal() + theme(axis.text = element_blank(), axis.title = element_blank(), axis.text.x = element_blank(),
			legend.position='right',  legend.direction='vertical', legend.text=element_text(size=7), legend.title=element_text(size=7))

		# combine and save
		assembled <- egg::ggarrange(hm, ds, nrow=1, ncol=2, widths=c(8, 0.3))
		ggsave(file.path(out_dir, paste0("Pearson_jaccard_heatmap.pdf")), assembled, width = 7, height = 6)
}


## MAIN
clusters_exclude <- c("BMMC22_ID2_hi_myeloid_prog")
datasets_include <- c("BMMC22", "PBMC9", "Islets", "K562", "Other")
res_id <- "granular_2M"

correlation_cols <- c("Pearson", "Pearson_log1p", "Spearman", "jaccard")
for (i in seq_along(correlation_cols)) {
	plot_single_correlation(correlation_cols[i], res_id, datasets_include, clusters_exclude)	
}
plot_double_correlation(res_id, datasets_include, clusters_exclude)