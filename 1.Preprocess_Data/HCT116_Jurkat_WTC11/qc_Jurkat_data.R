shhh <- suppressPackageStartupMessages # It's a library, so shhh!

shhh({
	library(tidyverse)
	library(Matrix)
	library(data.table)
	library(ggplot2)
	library(ggpubr)
	library(scales)
	})

# function to process each file
process_fragment_file <- function(file_path, cells, sample_id) {
	frag <- fread(file_path, sep = "\t", skip = "#", header = FALSE) %>% 
		setNames(c("chr", "start", "end", "barcode", "reads")) %>%
		mutate(barcode = paste0(sample_id, "_", barcode)) %>%
		filter(barcode %in% cells)

	print(head(frag))
	return(frag)
}

get_cells <- function(seu_path) {
		shhh({
		library(Seurat)
		library(Signac)
		})

	seurat_data <- get(load(seu_path))
	s <- subset(seurat_data, subset = final_annotation == "jurkat")
	cells <- rownames(s@meta.data)
	print(head(cells))
	message("Number of cells: ", length(cells))

	return(cells)
}

combine_fragment_files <- function(fragment_list, cells, out_dir) {
	# iterate through the files and combine results
	all_fragments <- lapply(names(fragment_list), function(sample) {
		process_fragment_file(fragment_list[[sample]], cells, sample)}) %>% 
		rbindlist() %>% as.data.frame()

	# sort
	all_fragments <- all_fragments[order(all_fragments$chr, as.numeric(all_fragments$start)), ]

	# write  output
	cat("Writing fragments to the output file...\n")
	temp_unzipped <- file.path(out_dir, "combined_atac_fragments.tsv")
	fwrite(all_fragments, temp_unzipped, sep = "\t", col.names = FALSE, quote = FALSE)
	remove(all_fragments)

	# zip and index
	bgzipped <- file.path(out_dir, "combined_atac_fragments.tsv.gz")
	Rsamtools::bgzip(temp_unzipped, bgzipped)  
	Rsamtools::indexTabix(bgzipped, format = "bed")

	system(paste("rm -f", temp_unzipped))
}

make_seurat_object <- function(seu_path, fragment_path, out_seu, out_meta) {
	shhh({
		library(Seurat)
		library(Signac)
		library(EnsDb.Hsapiens.v86) 
		library(BSgenome.Hsapiens.UCSC.hg38)
		})

	seurat_data <- get(load(seu_path))
	#print(str(seurat_data))

	s <- subset(seurat_data, subset = final_annotation == "jurkat")

	# annotate RNA
	s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
	s[["percent.rbio"]] <- PercentageFeatureSet(s, pattern = "^RPS|^RPL")

	# get annotations for TSS enrichment
	annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
	seqlevelsStyle(annotations) <- 'UCSC'
	genome(annotations) <- "hg38"

	# replace ATAC but save counts 
	atac_counts <- s@assays$ATAC@counts
	cells <- rownames(s@meta.data)

	s[["ATAC"]] <- NULL
	chrom_assay <- CreateChromatinAssay(counts = atac_counts, sep = c(":", "-"), genome = "hg38", fragments = fragment_path, min.cells = 10, annotation = annotations)
	s[["ATAC"]] <- chrom_assay
	DefaultAssay(s) <- "ATAC"

	# compute nuc signal and TSS enrichemnt
	s <- NucleosomeSignal(s)
	s <- TSSEnrichment(s)

	#str(s)
	saveRDS(s, out_seu)

	s_meta <- s@meta.data %>% 
		mutate(seurat_clusters = "Jurkat")
	fwrite(s_meta, out_meta, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
}

# make RNA qc plots qc plots
make_rna_qc_plots <- function(meta_path, cp, greys, rna_min, rna_max, gene_min, gene_max, pct_mt_max, pct_ribo_max, rna_qc_out) {
	# read in and format metadata
	meta <- read.csv(meta_path, sep = "\t", header = TRUE, row.names = 1)
	meta$cells <- rownames(meta)
	meta <- meta %>% 
		dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA) %>%
		mutate(nUMI_non_MT = nUMI*(100-percent.mt),
			percent.ribo = percent.rbio) %>%
		arrange(seurat_clusters) %>%
		mutate(seurat_clusters = as.character(seurat_clusters),
			seurat_clusters = replace_na(seurat_clusters, "NA"))
	meta$seurat_clusters <- factor(meta$seurat_clusters, levels = unique(meta$seurat_clusters), ordered = TRUE)

	# columns:  nUMI, nGene, percent.mt, percent.ribo, nCount_ATAC, nFeature_ATAC,
		# nucleosome_signal, nucleosome_percentile, TSS.enrichment, TSS.percentile, seurat_clusters,
		# prediction.score.max, predicted.id

	# plot first set of qc metrics
	pdf(rna_qc_out, width = 10, height = 12) 

	# number of cells per cluster
	cell_count <- table(meta$seurat_clusters) %>% as.data.frame()
	colnames(cell_count) <- (c("seurat_clusters", "n_cells"))

	p1 <- ggplot(cell_count, aes(x = seurat_clusters, y = n_cells)) +
		geom_bar(stat="identity", fill = greys[3]) + 
		geom_text(aes(label=n_cells), position=position_dodge(width=0.9), vjust=-0.25, size=3.5, angle = 45) +
		scale_fill_manual(values  = cp) +
		labs(title = "Cells per cluster", x = "Seurat cluster", y = "# cells") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position = "none")

	# distribution of UMIs
	p2 <- ggplot(meta, aes(x = nUMI)) + 
		geom_density(size=1, , fill = cp[10], color = cp[10], alpha = 0.5)  + 
		geom_vline(xintercept = c(rna_min, rna_max), color = greys[5], linetype = "dashed") +
		scale_x_log10(limits = c(100, 70000), n.breaks=10) + 
		labs(title = "UMI density", x = "# UMIs", y = "Cell density") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
			axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")
		
	# distribution of non-mito UMIs
	p3 <- ggplot(meta, aes(x = nUMI_non_MT)) + 
		geom_density(size=1, fill = cp[10], color = cp[10], alpha = 0.5)  + 
		scale_x_log10(limits = c(100, 70000), n.breaks=10) + 
		labs(title = "Non-mito UMI density", x = "# non-mito UMIs", y = "Cell density") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
			axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")

	# distribution of genes per cell
	p4 <- ggplot(meta, aes(x = nGene)) + 
		geom_density(size=1, fill = cp[10], color = cp[10], alpha = 0.5)  + 
		geom_vline(xintercept = c(gene_min, gene_max), color = greys[5], linetype = "dashed") +
		scale_x_log10(limits = c(100, 25000), n.breaks=10) + 
		labs(title = "# gene/cell distribution", x = "# genes per cell", y = "Cell density") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position = "none")
		
	# log-transformed genes per cell
	p5 <- ggplot(meta, aes(x=factor(seurat_clusters), y=log10(nGene))) +
		geom_boxplot(fill = cp[10], outlier.size = 1.5, outlier.shape = 16, outlier.color = greys[3]) +
		geom_hline(yintercept = log10(c(gene_min, gene_max)), color = greys[5], linetype = "dashed") +
		labs(title = "log10 (# genes/cell) per cluster", y = "log10 (# genes per cell)", x = "Seurat cluster") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position = "none")
	
	# gene counts versus UMI counts, colored by %mito
	p6 <- meta %>% arrange(percent.mt) %>%
		ggplot(aes(x = nUMI, y = nGene)) +
		geom_point(aes(color = percent.mt), shape = 16, size = 1) +
		geom_vline(xintercept = c(rna_min, rna_max), color = greys[5], linetype = "dashed") +
		geom_hline(yintercept = c(gene_min, gene_max), color = greys[5], linetype = "dashed") +
		scale_color_gradient(low = "#d3a9ce", high = "#430b4e") +
		labs(title = "Genes x UMIs, colored by % mito", x = "# UMIs per cell", y = "# genes per cell", color = "% mito") +
		stat_smooth(method = lm, color = greys[3]) +
		scale_x_log10(limits = c(100, 70000), n.breaks=10) +
		scale_y_log10(limits = c(100, 25000), n.breaks=10) +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
			axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), aspect.ratio = 1, legend.position = "right")

	# gene counts versus UMI counts, colored by %ribo
	p7 <- meta %>% arrange(percent.ribo) %>%
		ggplot(aes(x = nUMI, y = nGene)) +
		geom_point(aes(color = percent.ribo), shape = 16, size = 1) +
		geom_vline(xintercept = c(rna_min, rna_max), color = greys[5], linetype = "dashed") +
		geom_hline(yintercept = c(gene_min, gene_max), color = greys[5], linetype = "dashed") +
		scale_color_gradient(low = "#d3a9ce", high = "#430b4e") +
		labs(title = "Genes x UMIs, colored by % ribo", x = "# UMIs per cell", y = "# genes per cell", color = "% ribo") +
		stat_smooth(method = lm, color = greys[3]) +
		scale_x_log10(limits = c(100, 70000), n.breaks=10) +
		scale_y_log10(limits = c(100, 25000), n.breaks=10) +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
			axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), aspect.ratio = 1, legend.position = 'none')

	# nUMI versus % mito
	p7.5 <- meta %>% arrange(percent.ribo) %>%
		ggplot(aes(x=nUMI, y = percent.mt)) + 
		geom_point(aes(color = percent.ribo), shape = 16, size = 1) +
		geom_vline(xintercept = c(rna_min, rna_max), color = greys[5], linetype = "dashed") +
		labs(title = "UMIs x % mito, colored by % ribo", x = "# UMIs per cell", y = "% mitochondrial UMIs", color = "% ribo") +
		scale_color_gradient(low = "#d3a9ce", high = "#430b4e") +
		scale_x_log10(limits = c(100, 70000), n.breaks=10) +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
			axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), aspect.ratio = 1, legend.position = 'right')

	# distribution of mitochondrial  %
	p8 <- ggplot(meta, aes(x=percent.mt)) + 
		geom_density(size=1, fill = cp[10], color = cp[10], alpha = 0.5)  + 
		geom_vline(xintercept = c(pct_mt_max), color = greys[5], linetype = "dashed") +
		labs(title = "% mito density", x = "% mitochondrial UMIs", y = "Cell density") +
		scale_x_log10(n.breaks=8) +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position = "none")

	# distribution of ribosomal  %
	p9 <- ggplot(meta, aes(x=percent.ribo)) + 
		geom_density(size=1, fill = cp[10], color = cp[10], alpha = 0.5)  + 
		geom_vline(xintercept = c(pct_ribo_max), color = greys[5], linetype = "dashed") +
		labs(title = "% ribo density", x = "% ribosomal UMIs", y = "Cell density") +
		scale_x_log10(n.breaks=10) +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position = "none")

	# distribution of log10(genes/UMI) 
	p10 <- meta %>% mutate(log10_genes_per_UMI = log10(nGene/nUMI)) %>%
		ggplot(aes(x = log10_genes_per_UMI)) +
		geom_density(size=1, fill = cp[10], color = cp[10], alpha = 0.5)  + 
		labs(title = "log10 (genes per UMI)", x = "log10 (# genes / # UMIs)", y = "Cell density") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position = "none")

	# assemble and save
	figure <- ggarrange(ggarrange(p1, p2, p3, ncol = 3, widths = c(2, 1.5, 1.5)),
		ggarrange(p4, p5, ncol = 2, widths= c(1, 1.5)),
		ggarrange(p6, p7, p7.5, ncol = 3, widths = c(1.3, 1, 1.2)),
		ggarrange(p8, p9, p10, ncol = 3, widths=c(1, 1, 1)),
		nrow = 4, heights = c(1.2, 1, 1.3, 1))

	print(figure)
	dev.off()
}

# apply RNA filter
apply_rna_qc_filter <- function(meta_path, greys, rna_min, rna_max, gene_min, gene_max, pct_mt_max, pct_ribo_max, rna_filt_out) {
	# read in and format metadata
	meta <- read.csv(meta_path, sep = "\t", header = TRUE, row.names = 1)
	meta$cells <- rownames(meta)
	meta <- meta %>% 
		dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA) %>%
		mutate(nUMI_non_MT = nUMI*(100-percent.mt), percent.ribo = percent.rbio)

	meta_filt <- meta %>% dplyr::filter(!is.na(seurat_clusters),
		nUMI > rna_min, nUMI < rna_max,
		nGene > gene_min, nGene < gene_max,
		percent.mt < pct_mt_max, percent.ribo < pct_ribo_max)
		
	# number of cells per cluster after filtering
	cell_count <- table(meta_filt$seurat_clusters) %>% as.data.frame()
	colnames(cell_count) <- (c("seurat_clusters", "n_cells"))

	p1 <- ggplot(cell_count, aes(x = seurat_clusters, y = n_cells)) +
		geom_bar(stat="identity", fill = greys[3]) + 
		geom_text(aes(label=n_cells), position=position_dodge(width=0.9), vjust=-0.25, size=3.5, angle = 45) +
		scale_fill_manual(values  = cp) +
		labs(title = "Cells per cluster\nafter RNA QC filtering", x = "Seurat cluster", y = "# cells") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position = "none")
	ggsave(rna_filt_out, p1, width = 4, height = 3.5)
}

# make ATAC qc plots
make_atac_qc_plots_and_filter <- function(meta_path, greys, cp, atac_min, ataac_max, tss_enr_min, nuc_signal_max, atac_qc_out, atac_filt_out) {
	# read in and format metadata
	meta <- read.csv(meta_path, sep = "\t", header = TRUE, row.names = 1)
	meta$cells <- rownames(meta)
	meta <- meta %>% arrange(seurat_clusters) %>%
		mutate(seurat_clusters = as.character(seurat_clusters),
			seurat_clusters = replace_na(seurat_clusters, "NA"))
	meta$seurat_clusters <- factor(meta$seurat_clusters, levels = unique(meta$seurat_clusters), ordered = TRUE)
	
	pdf(atac_qc_out, width = 10, height = 7) 

	# number of cells per cluster
	cell_count <- table(meta$seurat_clusters) %>% as.data.frame()
	colnames(cell_count) <- (c("seurat_clusters", "n_cells"))

	# cells per cluster..
	p1 <- ggplot(cell_count, aes(x = seurat_clusters, y = n_cells)) +
		geom_bar(stat="identity", fill = greys[3]) + 
		geom_text(aes(label=n_cells), position=position_dodge(width=0.9), vjust=-0.25, size=3.5, angle = 45) +
		scale_fill_manual(values  = cp) +
		labs(title = "Cells per cluster", x = "Seurat cluster", y = "# cells") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position = "none")
	
	# distribution of ATAC fragments/cell
	p2 <- ggplot(meta, aes(x = nCount_ATAC)) + 
		geom_density(size=1, , fill = cp[7], color = cp[7], alpha = 0.5)  + 
		geom_vline(xintercept = c(atac_min, atac_max), color = greys[5], linetype = "dashed") +
		scale_x_continuous(n.breaks = 15) +
		#scale_x_log10(limits = c(100, 70000), n.breaks=10) + 
		labs(title = "Fragment density", x = "# fragments", y = "Cell density") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
			axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")

	# log10 ATAC fragments per cluster
	p3 <- ggplot(meta, aes(x=factor(seurat_clusters), y=log10(nCount_ATAC))) +
		geom_boxplot(fill = cp[7], outlier.size = 1.5, outlier.shape = 16, outlier.color = greys[3]) +
		geom_hline(yintercept = log10(c(atac_min, atac_max)), color = greys[5], linetype = "dashed") +
		labs(title = "log10 (# fragments/cell) per cluster", y = "log10 (# fragments per cell)", x = "Seurat cluster") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position = "none")

	# distribution of nucleosome signal 
	p4 <- ggplot(meta, aes(x = nucleosome_signal)) + 
		geom_density(size=1, , fill = cp[7], color = cp[7], alpha = 0.5)  + 
		geom_vline(xintercept = c(nuc_signal_max), color = greys[5], linetype = "dashed") +
		labs(title = "Nucleosome signal density", x = "Nucleosome signal", y = "Cell density") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position = "none")

	# distribution of TSS.enrichment per cluster
	p5 <- ggplot(meta, aes(x = seurat_clusters, y =TSS.enrichment)) + 
		geom_boxplot(fill = cp[7], outlier.size = 1.5, outlier.shape = 16, outlier.color = greys[3]) +
		geom_hline(yintercept = tss_enr_min, color = greys[5], linetype = "dashed") +
		ylim(c(2, 8)) +
		labs(title = "TSS enrichment per cluster", x = "Seurat cluster", y = "TSS enrichment") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position = "none")

	# # bin_breaks <- pretty(meta$TSS.enrichment, n = 7)  # Adjust `n` to control granularity
	# # num_bins <- length(bin_breaks) - 1  # Number of bins

	# contour_colors <- colorRampPalette(c("#ffffff", "#c5e5fb", "#002359"))(num_bins)

	# fragments x TSS enrichment
	num_bins <- 8
	contour_colors <- colorRampPalette(c("#ffffff", "#c5e5fb", "#002359"))(num_bins)

	# actually plot
	p6 <- ggplot(meta, aes(x = nCount_ATAC, y = TSS.enrichment)) + 
		geom_density_2d_filled(contour_var = "density", bins = num_bins) +
		geom_vline(xintercept = c(atac_min, atac_max), color = greys[5], linetype = "dashed") +
		geom_hline(yintercept = tss_enr_min, color = greys[5], linetype = "dashed") + 
		scale_fill_manual(values = contour_colors) +
		labs(title = "Fragments x TSS enrichment", x = "# fragments", y = "TSS enrichment", fill = "Density") +
		scale_x_log10(limits = c(100, 1e5), n.breaks=10) + 
		ylim(c(2, 8)) + 
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
			axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), aspect.ratio = 1,  legend.position = "right")

	# assemble and save
	figure <- ggarrange(ggarrange(p3, p5, ncol = 2, widths = c(1, 1)),
		ggarrange(p2, p4, p6, ncol = 3, widths= c(1, 1, 1.3)),
		nrow = 2, heights = c(1, 1))
	print(figure)
	dev.off()

	# plot number of cells after  ATAC filter
	meta_filt <- meta %>%
		dplyr::filter(seurat_clusters != "NA",
			nCount_ATAC > atac_min, nCount_ATAC < atac_max,
			nucleosome_signal < nuc_signal_max,
			TSS.enrichment > tss_enr_min) %>%
		mutate(seurat_clusters = as.numeric(seurat_clusters)) %>%
		dplyr::filter(seurat_clusters < 15)
		
	# number of cells per cluster after filtering
	cell_count <- table(meta_filt$seurat_clusters) %>% as.data.frame()
	colnames(cell_count) <- (c("seurat_clusters", "n_cells"))

	p1 <- ggplot(cell_count, aes(x = seurat_clusters, y = n_cells)) +
		geom_bar(stat="identity", fill = greys[3]) + 
		geom_text(aes(label=n_cells), position=position_dodge(width=0.9), vjust=-0.25, size=3.5, angle = 45) +
		scale_fill_manual(values  = cp) +
		labs(title = "Cells per cluster\nafter ATAC QC filtering", x = "Seurat cluster", y = "# cells") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position = "none")

	ggsave(atac_filt_out, p1, width = 4, height = 3.5)
}

# filter metadata by both RNA and ATAC metrics; save barcodes & stats and QC thresholds
apply_all_filters <- function(meta_path, greys, rna_min, rna_max, gene_min, gene_max, pct_mt_max, pct_ribo_max, atac_min, ataac_max, tss_enr_min, nuc_signal_max,
	cluster_stats_out, all_filt_out, cluster_barcodes_out, qc_thresholds_out) {
	meta <- read.csv(meta_path, sep = "\t", header = TRUE, row.names = 1)
	meta$cell <- rownames(meta)
	meta <- meta %>% 
		dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA) %>%
		mutate(nUMI_non_MT = nUMI*(100-percent.mt), percent.ribo = percent.rbio)

	# filter cells
	meta_filt <- meta %>% dplyr::filter(!is.na(seurat_clusters),
		nUMI > rna_min, nUMI < rna_max,
		nGene > gene_min, nGene < gene_max,
		percent.mt < pct_mt_max, percent.ribo < pct_ribo_max, 
		nCount_ATAC > atac_min, nCount_ATAC < atac_max,
		nucleosome_signal < nuc_signal_max,
		TSS.enrichment > tss_enr_min)

	# summarize clusters
	meta_summ <- meta_filt %>% group_by(seurat_clusters) %>%
		summarize(n_cells = n(), total_fragments = sum(nCount_ATAC), total_UMIs = sum(nUMI)) %>%
		mutate(mean_frag_per_cell = total_fragments / n_cells,
			mean_UMI_per_cell = total_UMIs / n_cells) %>%
		as_tibble()
	print(meta_summ, width = Inf)
	fwrite(meta_summ, cluster_stats_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

	# plot cells per cluster
	p1 <- ggplot(meta_summ, aes(x = factor(seurat_clusters), y = n_cells)) +
		geom_bar(stat="identity", fill = greys[4]) + 
		geom_text(aes(label=n_cells), position=position_dodge(width=0.9), vjust=-0.25, size=3.5, angle = 45) +
		scale_fill_manual(values  = cp) +
		labs(title = "Cells per cluster\nafter RNA+ATAC QC filtering", x = "Seurat cluster", y = "# cells") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position = "none")
	ggsave(all_filt_out, p1, width = 4, height = 3.5)

	# write barcode file
	meta_bc <- meta_filt %>% dplyr::select(cell, seurat_clusters)
	fwrite(meta_bc, cluster_barcodes_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

	# save qc thresholds for reference
	thresh <- data.frame(threshold = c("Minimum UMI per cell", "Maximum UMI per cell", "Minimum genes per cell", "Maximum genes per cell",
		"Maximum percent mitochondrial UMIs", "Maximum percent ribosomal UMIs", "Minimum fragments per cell", "Maximum fragments per cell",
		"Minimum TSS enrichment",  "Maximum nucleosomal signal"),
		value = c(rna_min, rna_max, gene_min, gene_max, pct_mt_max, pct_ribo_max, atac_min, ataac_max, tss_enr_min, nuc_signal_max))
	fwrite(thresh, qc_thresholds_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

}

filter_data_all_clusters <- function(seu_path, fragment_path, barcodes, out_dir) {
	shhh({
		library(Seurat)
		library(Signac)
	})

	seu <- readRDS(seu_path)
	cluster_bcs <- fread(barcodes) # cell, seurat_clusters
	clusters <- unique(cluster_bcs$seurat_clusters)
	fragments <- fread(fragment_path, header = FALSE) %>% 
		setNames(c("chr", "start", "end", "barcode", "counts"))
	message("Total fragments: ", nrow(fragments))

	#str(seu)

	for (c in seq_along(clusters)) {
		this_cluster <- clusters[c]

		message("Processing cluster ", this_cluster, "...")
		cluster_out <- file.path(out_dir, this_cluster); dir.create(cluster_out, showWarnings = FALSE)

		this_bc <- dplyr::filter(cluster_bcs, seurat_clusters == this_cluster) %>% pull(cell)
		message("# barcodes in metadata: ", length(this_bc))

		#seu_sub <- subset(seu, subset = final_annotation == this_cluster)
		save_cluster_fragment_file(fragments, this_bc, cluster_out)
		save_cluster_rna_matrix(seu, this_bc, cluster_out)
	}
}

save_cluster_fragment_file <- function(fragments, this_bc, out_dir) {
	all_fragments <- filter(fragments, barcode %in% this_bc)

	# sort
	all_fragments <- all_fragments[order(all_fragments[[1]], as.numeric(all_fragments[[2]])), ]
	num_frag <- nrow(all_fragments)
	message("# fragments in cluster: ", num_frag)

	# Write  output
	cat("Writing fragments to the output file...\n")
	temp_unzipped <- file.path(out_dir, "atac_fragments.tsv")
	fwrite(all_fragments, temp_unzipped, sep = "\t", col.names = FALSE, quote = FALSE)
	remove(all_fragments)

	bgzipped <- file.path(out_dir, "atac_fragments.tsv.gz")
	Rsamtools::bgzip(temp_unzipped, bgzipped)  
	Rsamtools::indexTabix(bgzipped, format = "bed")

	system(paste("rm -f", temp_unzipped))

}

save_cluster_rna_matrix <- function(seu_cluster, this_bc, out_dir) {
	rna_matrix <- seu_cluster@assays$RNA$counts
	colnames(rna_matrix) = dimnames(seu_cluster@assays$RNA$counts)[[2]] # cells
	rownames(rna_matrix) = dimnames(seu_cluster@assays$RNA$counts)[[1]] # genes

	rna_matrix <- rna_matrix[, this_bc]

	num_UMI <- sum(rna_matrix)
	message("# UMIs: ", num_UMI)

	rna_out <- file.path(out_dir, "rna_count_matrix.csv.gz")
	write.csv(as.data.frame(rna_matrix), gzfile(rna_out), quote = F)
}

### MAIN
## files
# fragment_path <- "/oak/stanford/groups/engreitz/Users/jray/231011-WTC11-V4-Multiome-ATAC/Cellranger_output/atac_fragments.tsv.gz"
# metadata_path <- "/oak/stanford/groups/engreitz/Users/ejagoda/cell_profiling/231011-WTC11-V4-e2g/cbc_w_seurat_clusters.txt" # skip line one and set colnames c("barcode", "cluster")
# cell_ranger_dir <- "/oak/stanford/groups/engreitz/Users/jray/231011-WTC11-V4-Multiome-ATAC/Cellranger_output"
# h5_cell_ranger_path <- "/oak/stanford/groups/engreitz/Users/jray/231011-WTC11-V4-Multiome-ATAC/Cellranger_output/filtered_feature_bc_matrix.h5"
# projected_metadata_path <- "/oak/stanford/groups/engreitz/Users/ejagoda/cell_profiling/pooled_wtc11_to_hpc_231011_w_scrublet_yes_filtered_w_jones_group_day_prediction.rds"

# for qc
out_dir <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0127_validation_cell_lines/Jurkat_processing"
# for saving per-cluster files
final_out_dir <-  "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0127_validation_cell_lines/scE2G_input"; dir.create(final_out_dir)

fragment_list <-c(e10l1 = "/oak/stanford/projects/igvf/Projects/250218_shared/e10l1_fragments.tsv.gz",
	e10l2 = "/oak/stanford/projects/igvf/Projects/250218_shared/e10l2_fragments.tsv.gz")
fragment_path <- file.path(out_dir, "combined_atac_fragments.tsv.gz")

seu_path <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0127_validation_cell_lines/data/seu_e10_final.Robj"
seu_out <- file.path(out_dir, "seu_Jurkat.rds")
meta_path <-  file.path(out_dir, "seu_Jurkat_metadata.tsv.gz")
rna_qc_out <- file.path(out_dir, "RNA_QC_plots_with_thresholds_v1.pdf")
rna_filt_out <- file.path(out_dir, "cells_per_cluster_after_RNA_QC.pdf")
atac_qc_out <- file.path(out_dir, "ATAC_QC_plots_with_thresholds_v1.pdf")
atac_filt_out <- file.path(out_dir, "cells_per_cluster_after_ATAC_QC.pdf")
cluster_stats_out <- file.path(out_dir, "filtered_cell_cluster_metrics.tsv")
all_filt_out <- file.path(out_dir, "cells_per_cluster_after_all_QC.pdf")
qc_thresholds_out <- file.path(out_dir, "qc_thresholds.tsv")
cluster_barcodes_out <- file.path(out_dir, "filtered_barcodes_with_clusters.tsv.gz")

## other params
colors <- c("#429130", "#0096a0", "#006eae", "#a64791", "#c5373d", "#e96a00", "#ca9b23")
cp <- c("#429130", "#2f9a71", "#159594", "#0096a0", "#0083ab",
	"#0075b3", "#006eae", "#5b5da3", "#8d4b9b", "#a64791", # 7 = blue, 10 = purple
	"#b03e67", "#c5373d", "#d8571f", "#e96a00", "#ca9b23")
greys <- c("#e5e5e9", "#c5cad7", "#96a0b3", "#6e788d", "#435369", "#1c2a43")

# thresholds (from wei-lin)
#min_cell <- 10
rna_min <- 10 ^ 3.6
rna_max <- Inf
pct_mt_max <- 30
pct_ribo_max <- 100
gene_min <- 0
gene_max <- Inf
atac_max <- Inf
atac_min <- 2e3
nuc_signal_max <- 1.5
tss_enr_min <- 3


### RUN

dir.create(out_dir, showWarnings = FALSE)
dir.create(final_out_dir, showWarnings = FALSE)

# cells <- get_cells(seu_path)
# combine_fragment_files(fragment_list, cells, out_dir)
#make_seurat_object(seu_path, fragment_path, seu_out, meta_path)  
# make_rna_qc_plots(meta_path, cp, greys, rna_min, rna_max, gene_min, gene_max, pct_mt_max, pct_ribo_max, rna_qc_out)
# apply_rna_qc_filter(meta_path, greys, rna_min, rna_max, gene_min, gene_max, pct_mt_max, pct_ribo_max, rna_filt_out)
# make_atac_qc_plots_and_filter(meta_path, greys, cp, atac_min, ataac_max, tss_enr_min, nuc_signal_max, atac_qc_out, atac_filt_out)
apply_all_filters(meta_path, greys, rna_min, rna_max, gene_min, gene_max, pct_mt_max, pct_ribo_max, atac_min, atac_max, tss_enr_min, nuc_signal_max,
	cluster_stats_out, all_filt_out, cluster_barcodes_out, qc_thresholds_out)

filter_data_all_clusters(seu_out, fragment_path, cluster_barcodes_out, final_out_dir) 