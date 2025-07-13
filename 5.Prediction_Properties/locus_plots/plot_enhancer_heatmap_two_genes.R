# attach required packages and functions
suppressPackageStartupMessages({
library(tidyverse)
library(data.table)
library(rtracklayer)
library(GenomicAlignments)
library(Gviz)
library(GenomicInteractions)
library(BiocParallel)})

make_enhancer_locus_plots <- function(row_to_run){
	##### INPUT 
	width = 6 # normal value: 8
	# params/files that are semi-permanent
	cluster_biosamples = TRUE # use gex-based clustering? if false, just sorts by n_enhancers
	project_dir = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0126_new_locus_plots"
	e2g_threshold <- 0.177 # for scE2G_Multiome
	param_scale_arc_width <- TRUE
	param_scale_arc_color <- TRUE
	n_enh_min <- 0 # minimum number of enhancers to include biosample
	pred_color <-  "#792374" # dark purple  #"#d3a9ce" # LIGHT purple 
	lightest_purple <- "#e9d3ea"
	darkest_purple <- "#430b4e"
	pred_out_color <-  "#792374" # dark purple  "#d3a9ce" #LIGHT purple
	atac_color <- "#006eae" # blue
	heatmap_colors <- c("#e9a0a5", "#730c0d") # red
	var_colors <- c("#f29742", "#e9c54e", "#bc9678", "#c5c1a5")
	sample_key_file = file.path(project_dir, "reference", "sample_key.tsv")
	gtf_path = "/oak/stanford/groups/engreitz/Users/sheth/Data/GENCODE/v43/gencode.v43.basic.annotation.gtf.gz"
	gene_ref_path = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_temp/scE2G/resources/genome_annotations/CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed"
	locus_key_file = file.path(project_dir, "reference", "locus_key.tsv")
	ordered_dim_file = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0214_new_global_properties/enhancer_similarity/granular_2M/ordered_biosample_names_Pearson_jaccard.rds"
	biosample_gex_clust_file = file.path(project_dir, "reference", "biosamples_tpm.tsv.gz")
	all_variants_file <- "/oak/stanford/groups/engreitz/Users/rosaxma/fine_mapped_UKBioBank/liftOver_hg38_GWAS_Traits/191010_UKBB_SUSIE/variant.list.all.tsv"
	cs_dir <- "/oak/stanford/groups/engreitz/Users/rosaxma/fine_mapped_UKBioBank/liftOver_hg38_GWAS_Traits/191010_UKBB_SUSIE/CredibleSetsFixed"


	## process inputs
	sample_key <- fread(sample_key_file)
	atac_files <- setNames(sample_key$ATAC_bw, sample_key$biosample) # biosample: atac_file
	biosample_names <- setNames(sample_key$biosample_name, sample_key$biosample) # biosample: biosample_name

	## get params from locus key
	locus_key <- fread(locus_key_file)

	plot_heatmap = locus_key$plot_heatmap[row_to_run] 
	plot_GWAS =  locus_key$plot_GWAS[row_to_run] 

	GENE_NAMES <- locus_key$gene[row_to_run] %>% strsplit(",") %>% unlist() %>% trimws()
	GENE_NAME <- paste(GENE_NAMES, collapse = "_")
	message("Locus: ", GENE_NAME)
	LOCUS <- gsub(",", "", locus_key$locus[row_to_run]) # remove commas from start+end
	LOCUS_CHR <- sub(":.*", "", LOCUS)
	positions <- sub(".*:", "", LOCUS)
	LOCUS_START <- as.numeric(sub("-.*", "", positions))
	LOCUS_END <- as.numeric(sub(".*-", "", positions))

	biosamples_expanded <- locus_key$biosamples_expanded[row_to_run] %>% strsplit(",") %>% unlist() %>% trimws()
	biosample_names_expanded <- biosample_names[biosamples_expanded] # get plotting names

	# read in prediction scores
	pred_score_files <- file.path(project_dir, GENE_NAMES, "thresholded_e2g_scores.rds")
	names(pred_score_files) <- GENE_NAMES
	pred_scores_list <- lapply(pred_score_files, FUN = readRDS)

	# define output files
	locus_for_files = paste0(LOCUS_CHR, "_", LOCUS_START, "_",  LOCUS_END)
	dir.create(file.path(project_dir, GENE_NAME))
	detailed_view_out = file.path(project_dir, GENE_NAME, paste0(locus_for_files, "_detailed_view.pdf"))
	heatmap_view_out = file.path(project_dir, GENE_NAME, paste0(locus_for_files, "_n", n_enh_min, "_heatmap_view.pdf"))
	combined_view_out = file.path(project_dir, GENE_NAME, paste0(locus_for_files, "_n", n_enh_min, "_combined_detailed_heatmap_view.pdf"))
	meta_plots_out = file.path(project_dir, GENE_NAME, paste0(locus_for_files, "_metadata_plots.pdf"))
	meta_table_out = file.path(project_dir, GENE_NAME, paste0(locus_for_files, "_metadata.tsv"))
	color_bar_out =  file.path(project_dir, GENE_NAME, paste0(locus_for_files, "_n", n_enh_min, "_color_bar.pdf"))
	credible_sets_out = file.path(project_dir, GENE_NAME, paste0(locus_for_files, "_credible_set_subset.tsv"))
	variants_out = file.path(project_dir, GENE_NAME, paste0(locus_for_files, "_variants_subset.tsv"))
	ovl_pred_out = file.path(project_dir, GENE_NAME, paste0(GENE_NAME, "_", locus_for_files, "_all_overlapping_predictions.tsv"))

	##### CREATE DETAILED VIEW
	message("Creating detailed view.")
	message("LOCUS_START: ", LOCUS_START, " LOCUS_END: ", LOCUS_END)

	## make genome coordinates track
	gtrack <- GenomeAxisTrack(scale = 0.85, label_pos = "left", lwd = 0.5, col = "#000000", fontcolor = "#000000")

		## make genes track
	# format genome annotations for gene track
	annot <- import(gtf_path, format = "gtf")
	gene_bed <- fread(gene_ref_path, sep= "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "Ensembl_ID", "gene_type"))

	# extract annotations for genes in universe and format as required
	genes <- annot %>% as.data.frame()
	genes <- genes %>% mutate(Ensembl_ID = sub("\\.\\d+$", "", gene_id)) %>% 
		filter(Ensembl_ID %in% gene_bed$Ensembl_ID, type == "exon") %>%
		select(chromosome = seqnames, start, end, width, strand, feature = gene_type, gene = gene_id,
				exon = exon_id, transcript = transcript_id, symbol = gene_name)

	# make track
	geneTrack = GeneRegionTrack(genes, genome = "hg38", name = "Genes",
								transcriptAnnotation = "symbol",
								collapseTranscripts = "meta",
								col = "#000000", fill = "#000000", col.line = "#C5C9D6",
								lwd = 0.5, cex = 0.3, fontface.group = 1, 
								chromosome = LOCUS_CHR, fontcolor = "#000000", fontcolor.group = "#000000",
								stacking = "squish", stackHeight = 0.75)

	## make interaction + bw tracks for expanded cell types
	# get GRanges of TSS to make interaction tracks
	gene_loci <- annot[annot$type == "gene"]
	make_tss_granges <- function(gene_id, gene_loci) {
		gene_tss <- resize(gene_loci[gene_loci$gene_name == gene_id], width = 1, fix = "start")
		return(gene_tss)
	}
	gene_tss <- lapply(GENE_NAMES, FUN = make_tss_granges, gene_loci = gene_loci)
	names(gene_tss) <- GENE_NAMES

	# functions to scale arc width and transparency based on score
	scale_arc_width <- function(scores, threshold, min_width = 0, max_width = 1, log10 = TRUE) {
		if (log10) {
			threshold <- log10(threshold)
			scores <- log10(scores)
		}

		scaled_scores <- (scores - threshold) / (1 - threshold)
		scaled_scores[scaled_scores < 0] <- 0  # Avoid negative values
		scaled_widths <- min_width + (max_width - min_width) * scaled_scores
		flipped <- 1 - scaled_widths
		return(flipped)

	}

	# Function to scale arc color based on score
	scale_arc_color <- function(scores, min_color, max_color, threshold, log10 = FALSE) {
		if (log10) {
			threshold <- log10(threshold)
			scores <- log10(scores)
		}

		# Scale scores to be between 0 and 1
		scaled_scores <- (scores - threshold) / (1 - threshold)
		scaled_scores[scaled_scores < 0] <- 0  # Clamp values below threshold to 0
		scaled_scores[scaled_scores > 1] <- 1  # Clamp values above 1 to 1

		# Convert hex colors to RGB (dividing by 255 to normalize to 0-1 range)
		min_rgb <- grDevices::col2rgb(min_color) / 255
		max_rgb <- grDevices::col2rgb(max_color) / 255
		# Replicate RGB values across scores for element-wise operations
		min_rgb <- matrix(min_rgb, nrow = 3, ncol = length(scores), byrow = FALSE)
		max_rgb <- matrix(max_rgb, nrow = 3, ncol = length(scores), byrow = FALSE)

		# Interpolate between min and max colors for each score
		interpolated_rgb <- min_rgb + (max_rgb - min_rgb) * matrix(scaled_scores, nrow = 3, ncol = length(scores), byrow = TRUE)

		# Convert interpolated RGB back to hex format
		interpolated_hex <- apply(interpolated_rgb, 2, function(rgb) {
			grDevices::rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 1)
		})
		
		return(interpolated_hex)
	}

	# function to process one biosample, returns c(pred_track, atac_track)
	make_expanded_tracks <- function(biosample_name, pairs, atac_file, locus_chr, locus_start, locus_end, atac_color){
		
		pred_track  <- InteractionTrack(pairs, chromosome = locus_chr, start = locus_start, end = locus_end, name = biosample_name)

		displayPars(pred_track) <- list(
			plot.outside = FALSE, 
			col.interactions = pairs$colors,
			col.anchors.line = pairs$colors,
			col.anchors.fill = pairs$colors,
			interaction.dimension = "width",
			interaction.measure = "fdr"
		)

		# atac track
		atac_track <- DataTrack(range = atac_file, genome = "hg38", name = "ATAC-seq",
			type = "polygon", fill.mountain = rep(atac_color, 2), lwd.mountain = 0, showAxis=FALSE)

		return(c(pred_track, atac_track))
	}

	make_interactions_across_genes_for_biosample <- function(pred_scores_list, biosample_id, gene_tss, GENE_NAMES,
		e2g_threshold, param_scale_arc_color, param_scale_arc_width) {
		pred_list <- list()
		tss_list <- list()

		for (i in seq_along(GENE_NAMES)){ # iterate across genes
			this_pred <- pred_scores_list[[i]][[biosample_id]]			
			pred_list <- c(pred_list, this_pred) 

			this_tss <- gene_tss[[i]]
			tss_list <- c(tss_list, this_tss[rep_len(1, length(this_pred))])
		}

		all_pred <- unlist(GRangesList(pred_list))
		all_tss <- unlist(GRangesList(tss_list))

		# sort and adjust color/width params
		all_pred <- all_pred[order(all_pred$score, decreasing = FALSE)]
		all_tss <- all_tss[order(all_pred$score, decreasing = FALSE)]

		if (param_scale_arc_color) {
			lightest_purple <- "#d3a9ce" # #b778b3  #e9d3ea #d3a9ce
			darkest_purple <- "#430b4e"
			scaled_colors <- scale_arc_color(all_pred$score, lightest_purple, darkest_purple, e2g_threshold)
		} else {
			scaled_colors <- pred_color	
		}
		
		if (param_scale_arc_width) {
			scaled_widths <- scale_arc_width(all_pred$score, e2g_threshold, 0, 2)
		} else {
			scaled_widths <- 1
		}

		pairs <- GenomicInteractions(all_pred, all_tss,
			fdr = scaled_widths, colors = scaled_colors)

		return(pairs)
	}

	# make the expanded tracks
	selected_biosample_tracks = list()
	selected_track_sizes = c()
	for (i in seq_along(biosamples_expanded)){
		biosample_id = biosamples_expanded[i]; print(biosample_id)
		this_pairs <- make_interactions_across_genes_for_biosample(pred_scores_list, biosample_id, gene_tss, GENE_NAMES,
			e2g_threshold, param_scale_arc_color, param_scale_arc_width)
		atac_file = atac_files[biosample_id]

		new_tracks = make_expanded_tracks(biosample_name = biosample_names_expanded[i], 
			pairs = this_pairs, atac_file = atac_file, 
			locus_chr = LOCUS_CHR, locus_start = LOCUS_START, locus_end = LOCUS_END,
			atac_color = atac_color)

		selected_biosample_tracks = c(selected_biosample_tracks, new_tracks)
		selected_track_sizes <- c(selected_track_sizes, 0.25, 0.2)
	}

	## combine "header" tracks for detailed view
	header_only_tracks <- c(gtrack, geneTrack)
	header_only_track_sizes <- c(0.15, 0.45) * 10

	header_tracks <- header_only_tracks
	header_track_sizes <- header_only_track_sizes

	## add two tracks for GWAS variants: one with dot plot of each variant x (PIP), one with credible sets (labelled by trait)
	if (plot_GWAS){
		message("Creating GWAS tracks.")
		traits <- locus_key$traits[row_to_run] %>% strsplit(",") %>% unlist() %>% trimws()
		select_var <-  locus_key$variants[row_to_run] %>% strsplit(",") %>% unlist() %>% trimws()

		# function to read in credible set file per trait
		read_gwas_file <- function(trait, this_chr, file_name, cs_dir){
		cs_file <- file.path(cs_dir, trait, file_name)
		cs <- fread(cs_file, sep="\t", header=TRUE) %>%
			dplyr::filter(chr == this_chr) %>%
			mutate(trait_CredibleSet = paste0(Disease, "_", CredibleSet))

		if (file_name == "variant.list.txt"){
			cs <- dplyr::filter(cs, pip > 0.1) %>%
				mutate(!!sym(trait) := pip)
			}
		return(cs)
		}

		# make track with individual variants
		var_df <- lapply(traits, read_gwas_file, this_chr = LOCUS_CHR, file_name = "variant.list.txt", cs_dir = cs_dir) %>%
			rbindlist(fill = TRUE) %>% as_tibble()

		var_gr <- makeGRangesFromDataFrame(var_df, seqnames.field = "chr", start.field = "start", end.field = "end",
			keep.extra.columns = TRUE, ignore.strand = TRUE, starts.in.df.are.0based = TRUE) %>%
			subsetByOverlaps(GRanges(seqnames = LOCUS_CHR, ranges = IRanges(LOCUS_START, LOCUS_END)))
		
		pip_df <- mcols(var_gr) %>% as_tibble() %>% dplyr::select(all_of(traits))

		var_track <- DataTrack(var_gr, type = c("p"), genome="hg38",
			data = pip_df, groups = traits, col=var_colors[1:length(traits)],
			cex = 0.5, ylim = c(0, 1), name = "GWAS variants\n(PIP)",
			fontcolor.legend = "#000000", fontface.legend = 1, fontsize.legend = 6)

		cs_df <- lapply(traits, read_gwas_file, this_chr = LOCUS_CHR, file_name = "all.cs.txt", cs_dir = cs_dir) %>%
			rbindlist() %>% as_tibble() %>%
			dplyr::filter(trait_CredibleSet %in% unique(var_df$trait_CredibleSet))
		cs_gr <- makeGRangesFromDataFrame(cs_df, keep.extra.columns = TRUE, ignore.strand = TRUE, starts.in.df.are.0based = TRUE) %>%
			subsetByOverlaps(GRanges(seqnames = LOCUS_CHR, ranges = IRanges(LOCUS_START, LOCUS_END)))

		# cs_track <- AnnotationTrack(cs_gr, genome = "hg38", shape = "box", id = mcols(cs_gr)$Disease, showFeatureId = TRUE,
		# 	fill = "#96a0b3", col = "#96a0b3",  name = "Credible sets", fontcolor.feature = "#000000",  cex = 0.3)
		cs_loc <- dplyr::select(cs_df, trait_CredibleSet, chr, cs_start = start, cs_end = end)
		var_loc <- dplyr::select(var_df, rsid, trait_CredibleSet, trait, chr, var_start = start, var_end = end) %>%
			mutate(strand = "*") %>%
			dplyr::filter(trait_CredibleSet %in% mcols(var_gr)$trait_CredibleSet)
		trait_col <- var_colors[1:length(traits)] %>% as.list() %>%  setNames(traits)

		#var_cs <- left_join(var_loc, cs_loc, by="trait_CredibleSet")
		cs_track <- AnnotationTrack(start = var_loc$var_start, end = var_loc$var_end, chromosome = LOCUS_CHR, #exon = var_loc$rsid,
			strand = var_loc$strand, group = var_loc$trait_CredibleSet, feature = var_loc$trait, #symbol = var_loc$trait,
			genome = "hg38",  shape = "box", fill = "#000000", col = "#000000", fontcolor = "#000000", cex = 0.2, lwd = 0.5,
			groupAnnotation = "feature", fontface.group = 1, min.width = 1,
			#!!!(trait_col), # need rlang
			showId = TRUE, #transcriptAnnotation = "symbol", 
			stackHeight = 0.75, name = "Credible sets")

		# add to detailed track list
		var_size <- ifelse(length(traits) == 1, 0.3, 0.6)
		header_tracks <- c(header_tracks, cs_track, var_track)
		header_track_sizes <- c(header_track_sizes, 0.3 * 10, var_size * 10)

		# make df with all predictions intersecting variant(s)
		if (length(select_var) > 0){
			all_pred_scores <- lapply(pred_scores_list, GRangesList) %>% lapply(unlist) %>% GRangesList() %>% unlist()
			ovl_pred_list <- list()
			for (i in seq_along(select_var)) {
				this_var <- select_var[i]
				select_var_gr <- var_gr[mcols(var_gr)$rsid == this_var]
				this_ovl <- subsetByOverlaps(all_pred_scores, select_var_gr, ignore.strand = TRUE) %>%
					as.data.frame() %>%
					mutate (rsid_overlap = this_var) %>%
					arrange(desc(score))
				ovl_pred_list[[length(ovl_pred_list) + 1]] <- this_ovl
			}
		ovl_pred <- rbindlist(ovl_pred_list) %>% as_tibble() %>% arrange(desc(score))
		fwrite(ovl_pred, ovl_pred_out, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

		}

		fwrite(as.data.frame(cs_gr) %>% arrange(start), credible_sets_out, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
		fwrite(as.data.frame(var_gr) %>% arrange(start), variants_out, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

		message("Added GWAS info.")

		ht = sum(header_track_sizes / 10) * 6 / 4 %>% round(2) # convert to inches 
		pdf(gwas_view_out, height = ht, width = width)
			plotTracks(trackList = header_tracks, chromosome = LOCUS_CHR, from = LOCUS_START, to = LOCUS_END,
				sizes = header_track_sizes,
				background.title = "transparent",
				col.axis = "#000000", fontcolor.title = "#000000", fontface.title = 1,
				rotation.title = 0, cex.axis = 0.3, cex.title = 0.5, margin = 30)
		dev.off()
		message("Plotted GWAS info")
		
	}


	# add expanded biosample predictions
	detailed_tracks <- c(header_only_tracks, selected_biosample_tracks)
	detailed_track_sizes <- c(header_only_track_sizes, selected_track_sizes * 10)
	message("Detailed track sizes: ", paste(detailed_track_sizes, collapse = ", "))

	## plot detailed view
	ht = sum(detailed_track_sizes / 10) * 5 / 4 %>% round(2) # convert to inches 
	pdf(detailed_view_out, height = ht, width = width)
		plotTracks(trackList = detailed_tracks, chromosome = LOCUS_CHR, from = LOCUS_START, to = LOCUS_END,
			sizes = detailed_track_sizes,
			background.title = "transparent",
			col.axis = "#000000", fontcolor.title = "#000000", fontface.title = 1,
			rotation.title = 0, cex.axis = 0.3, cex.title = 0.5, margin = 30)
	dev.off()
	message("Plotted detailed view.")

}

##### MAIN
rows_to_run = c(5)

for (i in seq_along(rows_to_run)){
	make_enhancer_locus_plots(rows_to_run[i])
}