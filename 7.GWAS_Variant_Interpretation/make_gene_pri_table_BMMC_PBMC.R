suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(ggdist)
})

plot_enrichment_distributions <- function(er_file, traits, cell_types, enr_threshold, out_dir) {
	# header: trait	nVariantsOverlappingEnhancers	biosample	nVariantsTotal	bpEnhancers	nCommonVariantsOverlappingEnhancers	nCommonVariantsTotal
	# recall	enrichment	method	biosampleGroup	traitGroup	SE_log_enr	CI_enr_low	CI_enr_high	p_enr	p_adjust_enr
	# recall_adjust	SE_recall	CI_recall_low	CI_recall_high

	er <- fread(er_file)  %>%
		mutate(group = "All traits & biosamples")
	selected <- dplyr::filter(er, biosample %in% cell_types, trait %in% traits) %>%
		mutate(group = "Selected traits & biosamples")
	sign <- dplyr::filter(er, p_adjust_enr < 0.05) %>%
		mutate(group = "Significantly-enriched")

	# plot distribution of enrichments
	to_plot <- rbind(er, selected, sign)
	group_counts <- group_by(to_plot, group) %>% summarize(n = n(), max_enr = max(enrichment)) %>% mutate(count_label = paste0("N = ", n))
	enr_lim = 21
	dist <- ggplot(to_plot, aes(x = enrichment, y = group)) +
		geom_vline(xintercept = c(1, enr_threshold), linetype = "dashed", color = "#c5cad7") +
		stat_eye(side = "top", shape = 16, point_size = 2, slab_linewidth = 0, slab_fill = "#96a0b3") +
		geom_text(data = group_counts, aes(x = enr_lim - 2, y = group, label = count_label), size = 2.5) +
		labs(y = "", x = "Enrichment (GWAS variants versus common variants)") +
		xlim(c(0, enr_lim)) +
		theme_classic() + theme( axis.text = element_text(size = 7), axis.title = element_text(size = 8))
	ggsave(file.path(out_dir, "enrichment_distributions.pdf"), dist, width = 5, height = 3)
}

get_enriched_pairs <- function(er_file, traits, cell_types, enr_threshold, out_dir) {
	er <- fread(er_file) %>%
		dplyr::filter(biosample %in% cell_types, trait %in% traits, p_adjust_enr < 0.05, enrichment > enr_threshold) %>%
		mutate(key = paste0(trait, "_", biosample))

	fwrite(er, file.path(out_dir, paste0("significantly_enriched_pairs_threshold", enr_threshold, ".tsv")), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

	return(unique(er$key))
}

get_distal_noncoding_distance_rank <- function(var, tss) {
	# genes within 2Mb
	lim = 2e6 #
	tss = dplyr::filter(tss, chr==var$chr[1], TSS>(min(var$start)-lim), TSS<(max(var$start)+lim)) %>% as_tibble()

	res_list = vector("list", nrow(tss))
	for (i in 1:length(res_list)){
		distToTSS = abs(tss$TSS[i]-var$start)		
		res_list[[i]] = data.frame(TargetGene=tss$gene[i], bestDNCDistanceToTSS=min(distToTSS))
	}

	res = rbindlist(res_list) %>% as_tibble %>%
		mutate(distanceRank = min_rank(bestDNCDistanceToTSS))

	return(res)
}

get_tpm_nlinks_per_cluster <- function(cell_types, pred_path_base, threshold) {
	cell_types_core = gsub("Islets_", "", cell_types, fixed=TRUE)

	ct_list = vector("list", length(cell_types))
	for (i in 1:length(ct_list)) {
		ct <- cell_types[i]
		pred_path <- ifelse(grepl("PBMC", ct, fixed = TRUE),
			paste0(pred_path_base, "PBMC"), paste0(pred_path_base, "BMMC"))
		file_path = file.path(pred_path, cell_types[i], "multiome_powerlaw_v2",  "encode_e2g_predictions.tsv.gz")
		pred =  fread(file_path, header=TRUE, sep="\t")
		colnames(pred)[colnames(pred)=="E2G.Score.qnorm"] = "pred_score"

		tpm = dplyr::select(pred, TargetGene, RNA_pseudobulkTPM) %>%
			distinct() %>%
			mutate(biosample=cell_types[i],
				RNA_pseudobulkTPM = round(RNA_pseudobulkTPM, 1))
		
		n_links = dplyr::filter(pred, pred_score>=threshold) %>%
			dplyr::select(TargetGene) %>% group_by(TargetGene) %>%
			summarize(nEnhancersPerGene = n())

		ct_list[[i]] = left_join(tpm, n_links, by="TargetGene")
	}

	ct_tpm = rbindlist(ct_list) %>% as_tibble
	summ = group_by(ct_tpm, TargetGene) %>%
		summarize(sumEnhancersPerGene = sum(nEnhancersPerGene, na.rm=TRUE))
	ct_tpm = left_join(ct_tpm, summ, by="TargetGene")

	return(ct_tpm)
}


get_results_per_credible_set <- function(intersect, var, tss, CredibleSet_this, trait_this) {
	intersect = dplyr::filter(intersect, CredibleSet==CredibleSet_this, trait==trait_this)
	var = dplyr::filter(var, CredibleSet==CredibleSet_this, trait==trait_this)
	distance_res = get_distal_noncoding_distance_rank(var, tss)

	res = intersect %>%
  		group_by(CredibleSet, trait, TargetGene) %>% # make one row per CS/trait/gene tuple
  		summarize(
			bestE2GVariantScore = max(predScore), # highest score value for the gene
			bestE2GCellType = biosample[which.max(predScore)], # biosample corresponding to the highest score
			bestE2GVariantID = rsid[which.max(predScore)], # variant with highest score
			nCellTypes = n_distinct(biosample),  # number of unique biosamples for the gene
			listCellTypes = paste(unique(biosample), collapse = ","), # comma-separated list of biosamples
			anyTPMOver1 = any(RNA_pseudobulkTPM>1), # T/F if any cell type has TPM>1
			listTPMs = paste(unique(RNA_pseudobulkTPM), collapse = ",")) %>% # comma-separated list of tpms
		ungroup() %>% group_by(CredibleSet, trait) %>% # now calculate rank for each group of CS/trait
		mutate(E2GRank = min_rank(desc(bestE2GVariantScore))) %>%
		ungroup() %>%
		left_join(distance_res, by="TargetGene") # add bestDNCDistanceToTSS and distanceRank 

	return(res)
}

# add pip of best E2G variant
add_pip <- function(df, all_var_file) {
	all_var = fread(all_var_file, header=FALSE, sep="\t") %>%
		setNames(c('chr', 'start', 'end', 'rsid', 'rsid2', 'ref', 'alt', 'trait', 'bestE2GVariantPIP')) %>%
		dplyr::select(rsid, trait, bestE2GVariantPIP)

	df = dplyr::left_join(df, all_var, by=c('bestE2GVariantID'='rsid', 'trait'='trait'))

	return(df)
}

merge_cs_annotations <- function(df, traits, cs_path) {
	# concat relevant all.cs.txt files
	all_cs_list = vector("list", length(traits))
	for (i in 1:length(all_cs_list)){
		file_path = file.path(cs_path, traits[i], "all.cs.txt")
		all_cs_list[[i]] = fread(file_path, header=TRUE, sep="\t") %>%
			dplyr::select(CredibleSet, LocusID, nSNP, AnyCoding, AnyPromoter, AnySpliceSite, Disease, BestSNP, BestSNPNearestGene) %>%
			dplyr::rename(trait=Disease, BestOverallSNP=BestSNP, BestSNPOverallNearestGene=BestSNPNearestGene)
	}
	all_cs = rbindlist(all_cs_list) %>% as_tibble

	# join with df
	df <- separate_wider_delim(df, CredibleSet, delim = "-", names = c("chr_start_hg19", "end_hg19", "cs_number"), cols_remove = FALSE) %>%
		mutate(LocusID = paste0(chr_start_hg19, "-", end_hg19)) %>%
		select(-c(chr_start_hg19, end_hg19, cs_number)) %>%
		left_join(all_cs, by=c("CredibleSet", "trait"))

	return(df)

}



# one row for each credible set-gene pair linked by scE2G
# columns: CredibleSet; trait; E2GRank; bestVariantE2GScore; bestE2GCellType; nCellTypes; listCellTypes; bestDistanceToTSS; distanceRank
make_gene_pri_table <- function(traits, cell_types, enriched_pairs, variant_file, pred_path_base, intersection_dir, tss_file, all_variants_file, cs_path, threshold, out_dir) {
	tss = fread(tss_file, header=FALSE, sep="\t") %>%
		setNames(c("chr", "start", "end", "gene", "score", "strand")) %>%
		mutate(TSS = (start+end)/2)

	# filter variants ot traits
	# chr	start	end	rsid	CredibleSet	trait
	var = fread(variant_file, header=TRUE, sep="\t") %>%
		dplyr::filter((trait %in% traits))

	# get tpm and num links data
	gene_metrics = get_tpm_nlinks_per_cluster(cell_types, pred_path_base, threshold)
	ct_gene_tpm = dplyr::select(gene_metrics, TargetGene, biosample, RNA_pseudobulkTPM)
	gene_nlinks = dplyr::select(gene_metrics, TargetGene, sumEnhancersPerGene) %>% distinct()

	# read in prediction/variant intersections, filter to traits, and merge across biosamples
	intersect_list = vector("list", length(cell_types))
	for (i in 1:length(cell_types)){
		#  /scratch/users/shethm/GWAS_E2G_benchmarking/results/2024_0711_sc_benchmarking_v1/scE2G_multiome/biosamples/Islets_PP/enhancerPredictions.thresholded.variantIntersection.tsv.gz
		# varChr	varStart	varEnd	rsid	CredibleSet	trait	predChr	predStart	predEnd	biosample	TargetGene	predScore
		file_name = file.path(intersection_dir, cell_types[i], "enhancerPredictions.thresholded.variantIntersection.tsv.gz")
		intersect_list[[i]] = fread(file_name, header=TRUE, sep="\t") %>%
			dplyr::filter(trait %in% traits) 
	}
	intersect = rbindlist(intersect_list) %>% as_tibble() %>%
		left_join(ct_gene_tpm, by=c("TargetGene", "biosample")) %>% # add RNA_pseudobulkTPM
		mutate(key = paste0(trait, "_", biosample)) %>%
		dplyr::filter(key %in% enriched_pairs)

	# iterate through unique credible set+trait pairs to get all E-G links and fill in...
	# E2G columns: bestVariantE2GScore; bestE2GCellType; nCellTypes; listCellTypes
	# distance columns: bestVariantDistanceToTSS, distanceRank [rank based on all variants]
	cs_vector = dplyr::select(intersect, CredibleSet, trait) %>% distinct()
	res_list = vector("list", nrow(cs_vector))
	for (i in 1:length(res_list)){
		res_list[[i]] = get_results_per_credible_set(intersect, var, tss, cs_vector$CredibleSet[i], cs_vector$trait[i])
	}
	res = rbindlist(res_list) %>% as_tibble() %>%
		add_pip(all_variants_file) %>%
		merge_cs_annotations(traits, cs_path) %>%
		mutate(isDistalNoncoding = (!AnyCoding & !AnyPromoter & !AnySpliceSite),
			isNearestGeneToBestSNP = (TargetGene==BestSNPOverallNearestGene)) %>%
		left_join(gene_nlinks, by="TargetGene") # add sumEnhancersPerGene

	# save
	dir.create(out_dir)
	fwrite(res, file.path(out_dir, "gene_prioritization_table_enr_filter.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, scipen=999)
	fwrite(intersect, file.path(out_dir, "extended_gene_prioritization_table_enr_filter.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, scipen=999)

	return(res)
}

make_filtered_table <- function(res, lim_cell_types, lim_dist, lim_e2g_rank, out_dir, dnc_vals = c(NA, TRUE)){
	filt <- dplyr::filter(res, nCellTypes <= lim_cell_types,
		bestDNCDistanceToTSS >= lim_dist, 
		E2GRank <= lim_e2g_rank,
		isDistalNoncoding %in% dnc_vals)

	fwrite(filt, file.path(out_dir, "prioritized_gene_prioritization_table_enr_filter.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, scipen=999)
}

### MAIN
## inputs to change
model_id = "scE2G_multiome"
GWAS_benchmark_id <- "2024_0911_hd_blood"
enr_threshold <- 5
intersection_dir = file.path("/scratch/users/shethm/GWAS_E2G_benchmarking/results", GWAS_benchmark_id, model_id, "biosamples")
pred_path_base = "/oak/stanford/groups/engreitz/Users/sheth/sc-E2G/results/2024_0826_"
er_file = file.path("/oak/stanford/groups/engreitz/Users/sheth/GWAS_benchmarking_working/GWAS_E2G_benchmarking/results",
	GWAS_benchmark_id, model_id, "variant_overlap", "enrichmentRecall.thresholded.traitByBiosample.tsv.gz")
out_dir ="/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0823_GWAS_variant_interpretation/BMMC22_PBMC9"
traits =  c("RBC", "WBC", "Mono", "Lym", "Eosino", "Baso", "Neutro", "Plt", "MCH", "MCHC", "MCV", "Hb", "HbA1c", "Ht", "Glucose")
cell_types_PBMC = c("PBMC9", "PBMC_10_cDC", "PBMC_11_CD14.Mono.1", "PBMC_12_CD14.Mono.2", "PBMC_13_CD16.Mono", "PBMC_17_B", "PBMC_20_CD4.N1", "PBMC_22_CD4.M", "PBMC_24_CD8.CM", "PBMC_25_NK")
cell_types_BMMC = c("BMMC22_B1_B", "BMMC22_CD4_pos_T_activated", "BMMC22_CD4_pos_T_naive", "BMMC22_CD8_pos_T", "BMMC22_CD8_pos_T_naive", "BMMC22_CD14_pos_Mono",
	"BMMC22_CD16_pos_Mono", "BMMC22_cDC2", "BMMC22_Erythroblast", "BMMC22_G_M_prog", "BMMC22_HSC", "BMMC22_ID2_hi_myeloid_prog", 
	"BMMC22_ILC", "BMMC22_Lymph_prog", "BMMC22_MK_E_prog", "BMMC22_Naive_CD20_pos_B", "BMMC22_NK", "BMMC22_Normoblast", "BMMC22_pDC", 
	"BMMC22_Plasma_cell", "BMMC22_Proerythroblast", "BMMC22_Transitional_B")
cell_types = c(cell_types_PBMC, cell_types_BMMC)

# inputs that will probably stay the same
tss_file = "/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-integrated/eQTLEnrichment/resources/genome_annotation/CollapsedGeneBounds.hg38.TSS500bp.bed"
variant_file = "/oak/stanford/groups/engreitz/Users/sheth/GWAS_benchmarking_working/GWAS_E2G_benchmarking/results/2024_0829_pancreas/variants/filteredGWASVariants.merged.sorted.tsv.gz"
all_variants_file = "/oak/stanford/groups/engreitz/Users/rosaxma/fine_mapped_UKBioBank/liftOver_hg38_GWAS_Traits/191010_UKBB_SUSIE/variant.list.all.tsv"
cs_path = "/oak/stanford/groups/engreitz/Users/rosaxma/fine_mapped_UKBioBank/liftOver_hg38_GWAS_Traits/191010_UKBB_SUSIE/CredibleSetsFixed/"
threshold = 0.164 # for scE2G_Multiome
# abc_anygene_file = "/oak/stanford/groups/engreitz/Users/sheth/GWAS_benchmarking_working/GWAS_E2G_benchmarking/resources/UKBiobank.ABCGene.anyabc.tsv"

## run
# plot_enrichment_distributions(er_file, traits, cell_types, 5, out_dir)
enriched_pairs <- get_enriched_pairs(er_file, traits, cell_types, enr_threshold, out_dir)
res <- make_gene_pri_table(traits, cell_types, enriched_pairs, variant_file, pred_path_base, intersection_dir,
	tss_file, all_variants_file, cs_path, threshold, out_dir)


