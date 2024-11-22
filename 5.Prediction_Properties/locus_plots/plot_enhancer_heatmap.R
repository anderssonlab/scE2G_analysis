# attach required packages and functions
suppressPackageStartupMessages({
library(tidyverse)
library(data.table)
library(rtracklayer)
library(rlang)
library(GenomicAlignments)
library(Gviz)
library(GenomicInteractions)
library(BiocParallel)})

make_enhancer_locus_plots <- function(row_to_run){
	##### INPUT 
	width = 6 # normal value: 8
	# params/files that are semi-permanent
	cluster_biosamples = TRUE # use gex-based clustering? if false, just sorts by n_enhancers
	project_dir = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0909_locus_plots"
	e2g_threshold <- 0.164 # for scE2G_Multiome
	n_enh_min <- 0 # minimum number of enhancers to include biosample
	pred_color <- "#d3a9ce" # LIGHT purple  "#792374" # dark purple 
	pred_out_color <- "#d3a9ce"
	atac_color <- "#006eae" # blue
	heatmap_colors <- c("#e9a0a5", "#730c0d") # red
	var_colors <- c("#f29742", "#e9c54e", "#bc9678", "#c5c1a5")
	sample_key_file = file.path(project_dir, "reference", "sample_key.tsv")
	gtf_path = "/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GENCODEv29/gencode.v29.annotation.gtf.gz"
	gene_ref_path = "/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/2023_1116.CollapsedGeneBounds.hg38.TSS500bp.bed"
	locus_key_file = file.path(project_dir, "reference", "locus_key.tsv")
	ordered_dim_file =  file.path(project_dir, "reference", "ordered_dimensions_granular.rds")
	biosample_gex_clust_file = file.path(project_dir, "reference", "biosamples_tpm.tsv.gz")
	all_variants_file <- "/oak/stanford/groups/engreitz/Users/rosaxma/fine_mapped_UKBioBank/liftOver_hg38_GWAS_Traits/191010_UKBB_SUSIE/variant.list.all.tsv"
	cs_dir <- "/oak/stanford/groups/engreitz/Users/rosaxma/fine_mapped_UKBioBank/liftOver_hg38_GWAS_Traits/191010_UKBB_SUSIE/CredibleSetsFixed"


	## process inputs
	sample_key <- fread(sample_key_file)
	atac_files <- setNames(sample_key$ATAC_bw, sample_key$biosample) # biosample: atac_file
	biosample_names <- setNames(sample_key$biosample_name, sample_key$biosample) # biosample: biosample_name

	## get params from locus key
	locus_key <- fread(locus_key_file, sep = "\t")

	plot_heatmap = locus_key$plot_heatmap[row_to_run] 
	plot_GWAS =  locus_key$plot_GWAS[row_to_run] 

	GENE_NAME <- locus_key$gene[row_to_run]; message("Locus: ", GENE_NAME)
	LOCUS <- gsub(",", "", locus_key$locus[row_to_run]) # remove commas from start+end
	LOCUS_CHR <- sub(":.*", "", LOCUS)
	positions <- sub(".*:", "", LOCUS)
	LOCUS_START <- as.numeric(sub("-.*", "", positions))
	LOCUS_END <- as.numeric(sub(".*-", "", positions))

	biosamples_expanded <- locus_key$biosamples_expanded[row_to_run] %>% strsplit(",") %>% unlist() %>% trimws()
	biosample_names_expanded <- biosample_names[biosamples_expanded] # get plotting names

	# read in prediction scores
	pred_score_file = file.path(project_dir, GENE_NAME, "thresholded_e2g_scores.rds")
	pred_scores <- readRDS(pred_score_file)

	# define output files
	locus_for_files = paste0(LOCUS_CHR, "_", LOCUS_START, "_",  LOCUS_END)
	detailed_view_out = file.path(project_dir, GENE_NAME, paste0(GENE_NAME, "_", locus_for_files, "_detailed_view.pdf"))
	heatmap_view_out = file.path(project_dir, GENE_NAME, paste0(GENE_NAME, "_", locus_for_files, "_n", n_enh_min, "_heatmap_view.pdf"))
	combined_view_out = file.path(project_dir, GENE_NAME, paste0(GENE_NAME, "_", locus_for_files, "_n", n_enh_min, "_combined_detailed_heatmap_view.pdf"))
	meta_plots_out = file.path(project_dir, GENE_NAME, paste0(GENE_NAME, "_", locus_for_files, "_metadata_plots.pdf"))
	meta_table_out = file.path(project_dir, GENE_NAME, paste0(GENE_NAME, "_", locus_for_files, "_metadata.tsv"))
	color_bar_out =  file.path(project_dir, GENE_NAME, paste0(GENE_NAME, "_", locus_for_files, "_n", n_enh_min, "_color_bar.pdf"))
	credible_sets_out = file.path(project_dir, GENE_NAME, paste0(GENE_NAME, "_", locus_for_files, "_credible_set_subset.tsv"))
	variants_out = file.path(project_dir, GENE_NAME, paste0(GENE_NAME, "_", locus_for_files, "_variants_subset.tsv"))
	ovl_pred_out = file.path(project_dir, GENE_NAME, paste0(GENE_NAME, "_", locus_for_files, "_all_overlapping_predictions.tsv"))

	##### CREATE DETAILED VIEW

	## make chromosome picture track
	ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = LOCUS_CHR)

	## make genome coordinates track
	gtrack <- GenomeAxisTrack(scale = 0.85, label_pos = "left", lwd = 0.5, col = "#000000", fontcolor = "#000000")

	## make genes track
	# format genome annotations for gene track
	annot <- import(gtf_path, format = "gtf")
	gene_bed = fread(gene_ref_path, header=FALSE, sep="\t") %>%
		setNames(c("chr", "start", "end", "name", "score", "strand"))

	# extract annotations on protein-coding genes and lincRNAs exons
	genes <- annot[annot$type == "exon" &
					annot$gene_type %in% c("protein_coding", "lincRNA") & 
					annot$transcript_type %in% c("protein_coding", "lincRNA")]

	# create data frame containing required gene annotation data
	genes <- genes %>% 
		as.data.frame() %>%
		select(chromosome = seqnames, start, end, width, strand, feature = gene_type, gene = gene_id,
				exon = exon_id, transcript = transcript_id, symbol = gene_name) %>%
			dplyr::filter(symbol %in% gene_bed$name) # filter to our gene universe

	# make track
	geneTrack = GeneRegionTrack(genes, genome = "hg38", name = "Genes",
								transcriptAnnotation = "symbol",
								collapseTranscripts = "meta",
								col = "#000000", fill = "#000000", # black
								lwd = 0.5, cex = 0.3, fontface.group = 1, 
								chromosome = LOCUS_CHR, fontcolor = "#000000",
								stacking = "squish", stackHeight = 0.75)

	## make interaction + bw tracks for expanded cell types
	# get GRanges of TSS to make interaction tracks
	gene_loci <- annot[annot$type == "gene"]
	gene_tss <- resize(gene_loci[gene_loci$gene_name == GENE_NAME], width = 1, fix = "start")

	# function to process one biosample, returns c(pred_track, atac_track)
	make_expanded_tracks <- function(biosample_name, pred_scores, atac_file, locus_chr, gene_tss, pred_color, pred_out_color, atac_color){
		# interaction track
		pairs <- GenomicInteractions(pred_scores,
			gene_tss[rep_len(1, length(pred_scores))],
			p.value = 1 - pred_scores$score - e2g_threshold)
		pred_track  <- InteractionTrack(pairs, chromosome = locus_chr, name = biosample_name)
		displayPars(pred_track) <- list(col.anchors.line = pred_color, col.anchors.fill = pred_color,
			plot.outside = TRUE, col.outside = pred_out_color,
			#col.interactions = alpha(pred_color, 0.7),
			col.interactions = pred_color,
			interaction.dimension = "width",
			interaction.measure = "p.value",
			interaction.dimension.transform = "log") 

		# atac track
		atac_track <- DataTrack(range = atac_file, genome = "hg38", name = "ATAC-seq",
							type = "polygon", fill.mountain = rep(atac_color, 2),
							lwd.mountain = 0, showAxis=FALSE,
							span = 1/6, evaluation = 100) # params for loess.smooth

		return(c(pred_track, atac_track))
	}

	# make the expanded tracks
	selected_biosample_tracks = list()
	for (i in seq_along(biosamples_expanded)){
		biosample_id = biosamples_expanded[i]; print(biosample_id)
		these_scores = pred_scores[[biosample_id]];
		atac_file = atac_files[biosample_id]

		new_tracks = make_expanded_tracks(biosample_name = biosample_names_expanded[i], 
			pred_scores = these_scores, atac_file = atac_file, gene_tss = gene_tss,
			locus_chr = LOCUS_CHR, pred_color = pred_color, pred_out_color = pred_out_color, atac_color = atac_color)

		selected_biosample_tracks = c(selected_biosample_tracks, new_tracks)
	}

	## combine "header" tracks for detailed view
	# header_tracks <- c(ideoTrack, gtrack, geneTrack)
	# header_track_sizes <- c(0.25, 0.4, 0.4)

	header_tracks <- c(gtrack, geneTrack)
	header_track_sizes <- c(0.15, 0.3)

	## add two tracks for GWAS variants: one with dot plot of each variant x -log10(PIP), one with credible sets (labelled by trait)
	if (plot_GWAS){
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
		print(var_df)

		var_gr <- makeGRangesFromDataFrame(var_df, seqnames.field = "chr", keep.extra.columns = TRUE,
			ignore.strand = TRUE, starts.in.df.are.0based = TRUE) %>%
			subsetByOverlaps(GRanges(seqnames = LOCUS_CHR, ranges = IRanges(LOCUS_START, LOCUS_END)))

		pip_df <- mcols(var_gr) %>% as_tibble() %>% dplyr::select(all_of(traits))
		var_track <- DataTrack(var_gr, type = c("p"), genome="hg38",
			data = pip_df, groups = traits, col=var_colors[1:length(traits)],
			cex = 0.5, ylim = c(0, 1), name = "GWAS variants\n(PIP)",
			fontcolor.legend = "#96a0b3", fontface.legend = 1, fontsize.legend = 6)

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
			genome = "hg38",  shape = "box", fill = "#6e788d", col = "#6e788d", fontcolor = "#6e788d", cex = 0.2, lwd = 0.5,
			groupAnnotation = "feature", fontface.group = 1, min.width = 1,
			#!!!(trait_col), # need rlang
			showId = TRUE, #transcriptAnnotation = "symbol", 
			stackHeight = 0.75, name = "Credible sets")

		# add to detailed track list
		var_size <- ifelse(length(traits) == 1, 0.3, 0.6)
		header_tracks <- c(header_tracks, cs_track, var_track)
		header_track_sizes <- c(header_track_sizes, 0.3, var_size)

		# make df with all predictions intersecting variant(s)
		if (length(select_var) > 0){
			all_pred_scores <-  unlist(GRangesList(pred_scores))
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
	}

	# add expanded biosample predictions
	detailed_tracks <- c(header_tracks, selected_biosample_tracks)
	detailed_track_sizes <- c(header_track_sizes, rep(c(0.3, 0.2), times = length(biosamples_expanded)))

	## plot detailed view
	ht = sum(detailed_track_sizes) * 5 / 4 # convert to inches 
	pdf(detailed_view_out, height = ht, width = width)
		plotTracks(trackList = detailed_tracks, chromosome = LOCUS_CHR, from = LOCUS_START, to = LOCUS_END,
			sizes = detailed_track_sizes,
			background.title = "transparent",
			col.axis = "#000000", fontcolor.title = "#000000", fontface.title = 1,
			rotation.title = 0, cex.axis = 0.3, cex.title = 0.5, margin = 30)
	dev.off()
	message("Plotted detailed view.")

	if (!plot_heatmap) {quit()}

	##### CREATE HEATMAP ACROSS ALL BIOSAMPLES
	message("Making heatmap...")
	## filter + format predictions
	# function to subset enhancers to selected region 
	get_links <- function(scores, chr, start, end) {
	region <- GRanges(seqnames = chr, ranges = IRanges(start, end))
	subsetByOverlaps(scores, ranges = region)
	}

	# filter scores for enhancers with selected locus region
	pred_scores <- lapply(pred_scores, FUN = get_links, chr = LOCUS_CHR, start = LOCUS_START, end = LOCUS_END)
	message("Read in all scores.")

	# filter by number of enhancers
	n_enh <- vapply(pred_scores, FUN = length, FUN.VALUE = integer(1)) # calculate n enh
	n_enh <- n_enh[n_enh >= n_enh_min] # filter
	pred_scores <- pred_scores[names(n_enh)]

	# sort biosamples
	if (cluster_biosamples){
		# sort based on gene expression clustering
		ordered_dim <- readRDS(ordered_dim_file)
		biosample_names_order <- ordered_dim[["biosamples_cor"]] # names not actual biosample!
		biosample_order <- sample_key$biosample[match(biosample_names_order, sample_key$biosample_name)] # get corresponding biosamples
		biosample_order <- intersect(biosample_order, names(pred_scores)) # filter to available biosamples
	} else {
		# sort by descending number of enhancers
		biosample_order <- names(n_enh)[rev(order(n_enh))] # sort descending
	}
	pred_scores <- pred_scores[biosample_order] # extract scores for these samples
	message("Sorted biosamples.")

	# function to combine GRanges for each biosample into single object with  metadata column per sample (h/t ChatGPT)
	combine_granges <- function(gr_list) {
	# find all unique ranges across all samples
	sample_names = names(gr_list)
	gr_list = GRangesList(gr_list)
	all_ranges <- reduce(unlist(gr_list))
	combined_gr <- all_ranges
	
	# iterate through samples
	for (i in seq_along(gr_list)) {
		sample_gr <- gr_list[[i]]
		overlaps <- findOverlaps(combined_gr, sample_gr) # get overlaps with sample
		
		score_column <- rep(NA, length(combined_gr)) # initialize score vector with NAs
		score_column[queryHits(overlaps)] <- sample_gr$score[subjectHits(overlaps)]
		
		mcols(combined_gr)[[sample_names[i]]] <- score_column # add this sample's scores as new metadata column
	}
	return(combined_gr)
	}

	# combine scores for all biosamples
	combined_scores <- combine_granges(pred_scores)
	message("Merged GRanges.")

	## create heatmap data track 
	e2g_tracks_combined <- DataTrack(combined_scores, type = c("heatmap"), showSampleNames = TRUE, 
				name = " ", separator = 1, genome = "hg38", ylim = c(e2g_threshold, 1), 
				gradient = heatmap_colors, showColorBar = FALSE, showAxis = FALSE,  na.rm=TRUE,
			cex.sampleNames = 0.3, col.sampleNames="#000000")
	message("Made main heatmap track.")

	## plotting
	# track sizes for 1) combined detailed + heatmap view; 2) just heatmap view
	heatmap_size = length(pred_scores) * 0.04
	combined_tracks <- c(detailed_tracks, e2g_tracks_combined)
	combined_track_sizes <- c(detailed_track_sizes, heatmap_size)
	heatmap_tracks <- c(header_tracks, e2g_tracks_combined)
	heatmap_track_sizes <- c(header_track_sizes, heatmap_size)

	# plot heatmap only
	ht = sum(heatmap_track_sizes) * 5 / 4
	pdf(heatmap_view_out, height = ht, width = width)
		plotTracks(trackList = heatmap_tracks, chromosome = LOCUS_CHR, from = LOCUS_START, to = LOCUS_END,
				sizes = heatmap_track_sizes,
				background.title = "transparent",
				col.axis = "#000000", fontcolor.title = "#000000", fontface.title = 1,
				rotation.title = 0, cex.axis = 0.3, cex.title = 0.5, margin = 30)
	dev.off()
	message("Plotted heatmap-only view!")

	# plot combined detailed + heatmap
	ht = sum(combined_track_sizes) * 5 / 4
	pdf(combined_view_out, height = ht, width = width)
		plotTracks(trackList = combined_tracks, chromosome = LOCUS_CHR, from = LOCUS_START, to = LOCUS_END,
				sizes = combined_track_sizes,
				background.title = "transparent",
				col.axis = "#000000", fontcolor.title = "#000000", fontface.title = 1,
				rotation.title = 0, cex.axis = 0.3, cex.title = 0.5, margin = 30)
	dev.off()
	message("Plotted combined detailed + heatmap view!")

	## plot heatmap with score color bar for illustrator
	colorbar = DataTrack(combined_scores, type = c("heatmap"), showSampleNames = FALSE, 
				name = "heatmap", separator = 1, genome = "hg38",
				gradient = heatmap_colors, showColorBar = TRUE, showAxis = TRUE, na.rm = TRUE)
	colorbar_tracks = c(gtrack, colorbar)
	colorbar_sizes = c(0.4, heatmap_size)

	pdf(color_bar_out, height = 8, width = width)
		plotTracks(colorbar_tracks, chromosome = LOCUS_CHR, from = LOCUS_START, to = LOCUS_END,
			sizes = colorbar_sizes, col.axis="#000000", 
			background.title = "transparent")
	dev.off()
	message("Plotted color bar.")

	##### RANDOM EXTRA PLOTS

	## some params
	greys <- c("#e5e5e9", "#1c2a43") 
	tpm_lims <- c(0, 100)

	## get TPM for gene per biosample
	tpm_df <- fread(biosample_gex_clust_file) %>%
		setNames(c("TargetGene", "biosample_name", "TPM")) %>%
		dplyr::filter(TargetGene == GENE_NAME) %>%
		dplyr::left_join(sample_key, by="biosample_name") %>%
		mutate(TPM = pmax(TPM, tpm_lims[1]), TPM = pmin(TPM, tpm_lims[2]))

	## get enhancer info and add TPM
	mdata <- mcols(combined_scores) %>% as_tibble() %>% # columns: biosamples; rows: enhancers
		pivot_longer(cols = everything(), names_to = "biosample", values_to = "score") %>%
		dplyr::filter(!is.na(score)) %>%
		group_by(biosample) %>%
		summarize(n_enh = n(),
			score_sum = sum(score),
			score_mean = mean(score),
			score_max = max(score)) %>%
		left_join(tpm_df, by="biosample")

	fwrite(mdata, meta_table_out, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
	message("Made & saved metadata df.")

	## format for plotting
	label_key <- c(
		"n_enh" = "# enhancers",
		"score_sum" = "Sum of scE2G scores",
		"score_mean" = "Mean scE2G score",
		"score_max" = "Maximum scE2G score")

	df_plot <- mdata %>%
		dplyr::select(TPM, biosample, biosample_name, all_of(names(label_key))) %>%
		pivot_longer(cols = -c(TPM, biosample, biosample_name), names_to = "metric", values_to = "value")
	df_exp <- dplyr::filter(df_plot, biosample %in% biosamples_expanded)

	## make scatterplot
	p <- ggplot(df_plot, aes(x = TPM, y = value)) +
		geom_point(alpha = 0.7, color = greys[2]) +
		geom_point(data = df_exp, aes(x = TPM, y = value), alpha = 0.8, size = 3, color = pred_color) +
		facet_wrap(~ metric, scales = "free_y", nrow = 1, labeller = as_labeller(label_key), strip.position = "left") + 
		labs(y = NULL) +
		#geom_text_repel(data = df_exp,  aes(label = biosample_name),  size = 3, color = "#000000", max.overlaps = 10) +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), 
				legend.position = "None", strip.background = element_blank(), strip.placement = "outside")

	ggsave(meta_plots_out, p, width=8, height=2)
	message("Plotted metadata!")

}



##### MAIN
rows_to_run = c(14,15)

for (i in seq_along(rows_to_run)){
	make_enhancer_locus_plots(rows_to_run[i])
}