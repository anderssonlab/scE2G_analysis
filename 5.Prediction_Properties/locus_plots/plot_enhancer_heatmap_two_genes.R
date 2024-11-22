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
	project_dir = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0909_locus_plots"
	e2g_threshold <- 0.164 # for scE2G_Multiome
	n_enh_min <- 0 # minimum number of enhancers to include biosample
	pred_color <- "#d3a9ce" # purple
	pred_out_color <- "#d3a9ce"
	atac_color <- "#006eae" # blue
	heatmap_colors <- c("#e9a0a5", "#730c0d") # reds
	var_colors <- c("#f29742", "#e9c54e", "#bc9678", "#c5c1a5")
	sample_key_file = file.path(project_dir, "reference", "sample_key.tsv")
	gtf_path = "/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GENCODEv29/gencode.v29.annotation.gtf.gz"
	gene_ref_path = "/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/2023_1116.CollapsedGeneBounds.hg38.TSS500bp.bed"
	locus_key_file = file.path(project_dir, "reference", "locus_key.tsv")
	ordered_dim_file =  file.path(project_dir, "reference", "ordered_dimensions.rds")
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
	# get GRanges of TSS to make interaction tracks
	gene_loci <- annot[annot$type == "gene"]
	make_tss_granges <- function(gene_id, gene_loci) {
		gene_tss <- resize(gene_loci[gene_loci$gene_name == gene_id], width = 1, fix = "start")
		return(gene_tss)
	}
	gene_tss <- lapply(GENE_NAMES, FUN = make_tss_granges, gene_loci = gene_loci)
	names(gene_tss) <- GENE_NAMES

	# function to process one biosample, returns c(pred_track, atac_track)
	make_expanded_tracks <- function(biosample_name, pairs, atac_file, locus_chr, pred_color, pred_out_color, atac_color){
		pred_track  <- InteractionTrack(pairs, chromosome = locus_chr, name = biosample_name)
		displayPars(pred_track) <- list(col.interactions = pred_color, col.anchors.line = pred_color, col.anchors.fill = pred_color,
			plot.outside = TRUE, col.outside = pred_out_color,
			interaction.measure = "p.value",
			#interaction.dimension.transform = "log",
			interaction.dimension = "width") # other option is width

		# atac track
		atac_track <- DataTrack(range = atac_file, genome = "hg38", name = "ATAC-seq",
							type = "polygon", fill.mountain = rep(atac_color, 2),
							lwd.mountain = 0, showAxis=FALSE)

		return(c(pred_track, atac_track))
	}

	make_interactions_across_genes_for_biosample <- function(pred_scores_list, biosample_id, gene_tss, GENE_NAMES, e2g_threshold) {
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

		pairs <- GenomicInteractions(all_pred, all_tss, p.value = 1 - all_pred$score)

		return(pairs)
	}

	# make the expanded tracks
	selected_biosample_tracks = list()
	for (i in seq_along(biosamples_expanded)){
		biosample_id = biosamples_expanded[i]; print(biosample_id)
		this_pairs <- make_interactions_across_genes_for_biosample(pred_scores_list, biosample_id, gene_tss, GENE_NAMES, e2g_threshold)
		atac_file = atac_files[biosample_id]

		new_tracks = make_expanded_tracks(biosample_name = biosample_names_expanded[i], 
			pairs = this_pairs, atac_file = atac_file, 
			locus_chr = LOCUS_CHR, pred_color = pred_color, pred_out_color = pred_out_color, atac_color = atac_color)

		selected_biosample_tracks = c(selected_biosample_tracks, new_tracks)
	}

	## combine "header" tracks for detailed view
	header_tracks <- c(gtrack, geneTrack)
	header_track_sizes <- c(0.15, 0.3)

		## add two tracks for GWAS variants: one with dot plot of each variant x (PIP), one with credible sets (labelled by trait)
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
			data = pip_df, groups = traits,
			col=var_colors[1:length(traits)], cex = 0.5, ylim = c(0, 1), name = "GWAS variants\n(PIP)")

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
			groupAnnotation = "feature", fontface.group = 1, #just.group = "above",
			#!!!(trait_col), # need rlang
			showId = TRUE, #transcriptAnnotation = "symbol", 
			stackHeight = 0.75, name = "Credible sets")

		# add to detailed track list
		header_tracks <- c(header_tracks, cs_track, var_track)
		header_track_sizes <- c(header_track_sizes, 0.4, 0.3)

		# make df with all predictions intersecting variant(s)
		if (length(select_var) > 0){
			all_pred_scores <- lapply(pred_scores_list, GRangesList) %>% lapply(unlist) %>% GRangesList() %>% unlist()
 			#all_pred_scores <-  unlist(GRangesList(pred_scores))
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
	detailed_track_sizes <- c(header_track_sizes, rep(c(0.25, 0.15), times = length(biosamples_expanded)))

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

}

##### MAIN
rows_to_run = c(16)

for (i in seq_along(rows_to_run)){
	make_enhancer_locus_plots(rows_to_run[i])
}