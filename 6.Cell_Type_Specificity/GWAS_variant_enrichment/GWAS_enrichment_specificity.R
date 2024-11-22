suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(ggdist)
  library(cowplot)
})

plot_enrichment_across_biosamples <- function(this_trait, er, out_dir) {
	filt <- dplyr::filter(er, trait == this_trait) %>%
		arrange(enrichment)
	filt$biosample <- factor(filt$biosample, levels = filt$biosample, ordered = TRUE)

	bar <- ggplot(filt, aes(x = enrichment, y = biosample)) + 
		geom_col(aes(fill = significant_enr)) +
		geom_errorbar(aes(xmin=CI_enr_low, xmax=CI_enr_high), width=0.2, linewidth = 0.25) +
		scale_fill_manual(values = c(not_significant="#b778b3", significant="#792374")) +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position = "none")
	
	ggsave(file.path(out_dir, paste0(this_trait, "_enrichment_all.pdf")), bar, height = 7, width = 3)
}

compare_enrichment_across_biosamples <- function(this_trait, top_n, er, out_dir) {
	filt <- dplyr::filter(er, trait == this_trait) %>%
		arrange((desc(enrichment)))
	top_biosamples <- filt$biosample[1:top_n]
	others <- filt$biosample[(top_n+1):(nrow(filt))]
	groups <- expand.grid(top_biosample = top_biosamples, other_biosample = others)

	# iterate through rows
	res_list <- vector("list", nrow(groups))
  	for (i in 1:nrow(groups)) {
		group1 <- filt %>% dplyr::filter(biosample == groups$top_biosample[i])
		group2 <- filt %>% dplyr::filter(biosample == groups$other_biosample[i])

		d = log(group1$enrichment[1]) - log(group2$enrichment[1]) # difference in log RRs
		SE_d = sqrt(group1$SE_log_enr[1]**2 + group2$SE_log_enr[1]**2)
		z = d/SE_d
		p = pnorm(-(abs(z))) # ONE-sided p-value
		 
		 res_list[[i]] <- as.data.frame(list(group1 = groups$top_biosample[i], group2 = groups$other_biosample[i],
		 	group1_enr= group1$enrichment, group2_enr=group2$enrichment, delta_enr=group1$enrichment-group2$enrichment,
			p_delta_enr = p))
  	}
	
	res <- rbindlist(res_list) %>% as_tibble() %>%
		mutate(p_delta_enr_adjust = p.adjust(p_delta_enr, method = "bonferroni"),
			significant = p_delta_enr_adjust < 0.05) %>%
		distinct() %>%
		arrange(desc(group2_enr))
	fwrite(res, file.path(out_dir, paste0(this_trait, "_enrichment_significance_top", top_n, ".tsv")), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	
	res <- res
}

plot_trait_summary <- function(this_trait, top_n, first_sign_n, er, out_dir) {
	er <- dplyr::filter(er, !(grepl("K562", biosample_name, fixed = TRUE)), !(grepl("GM12878", biosample_name, fixed = TRUE)), !is.na(biosample_name))
	res <- compare_enrichment_across_biosamples(this_trait, top_n, er, out_dir)
	
	filt <- dplyr::filter(er, trait == this_trait) %>%
		arrange((desc(enrichment))) %>%
		mutate(enr_rank = row_number())
		
	res_summ <- dplyr::filter(res, group2 %in% filt$biosample[first_sign_n:nrow(filt)]) %>%
		group_by(group1) %>%
		summarize(max_p_adjust = max(p_delta_enr_adjust),
			significant_delta = max_p_adjust < 0.05,
			stars = case_when(max_p_adjust < 0.001 ~ "***", max_p_adjust < 0.01 ~ "**", max_p_adjust < 0.05 ~ "*", max_p_adjust < 1 ~ ""))
		
	filt <- filt %>% left_join(res_summ, by = c("biosample"="group1")) %>%
		mutate(stars = replace_na(stars, ""), 
			stars_loc = CI_enr_high + 1,
			color = case_when(enr_rank <= top_n ~ "#430b4e", enr_rank < first_sign_n ~ "#a64791", enr_rank >= first_sign_n ~ "#96a0b3")) 	

	filt_top <- dplyr::filter(filt, enr_rank < first_sign_n)
	print(filt_top$stars)
	print(filt_top$stars_loc)

	#print(filt$biosample_name[first_sign_n:nrow(filt)])
	filt_other <- dplyr::filter(filt, enr_rank >= first_sign_n) %>%
		# mutate(biosample_name = case_when(grepl("BMMC", biosample_name) ~ "Other BMMC clusters",
		# 		grepl("PBMC", biosample_name) ~ "Other PBMC clusters",
		# 		grepl ("Pancreatic", biosample_name) ~ "Other islet clusters"))
		mutate(biosample_name = "Other biosample")
	filt_summ <- group_by(filt_other, biosample_name) %>% summarize(avg_enr = mean(enrichment)) %>% arrange(avg_enr)
	
	
	set_levels <- (c(filt_summ$biosample_name, rev(filt_top$biosample_name)))
		
	filt_top$biosample_name <- factor(filt_top$biosample_name, levels = set_levels, ordered = TRUE)
	filt_other$biosample_name <- factor(filt_other$biosample_name, levels = set_levels, ordered = TRUE)

	x_label <- paste0("Enrichment\nGWAS variants versus common variants\nN = ", filt$nVariantsTotal[1], " ", this_trait, " variants") 
	enr_max <- max(filt_top$CI_enr_high)

	top <- ggplot(filt_top, aes(x = enrichment, y = biosample_name)) + 
		geom_point(aes(fill = color, color = color), shape = 16, size = 3.5) +
		geom_errorbar(aes(xmin=CI_enr_low, xmax=CI_enr_high, color = color), width=0.2, linewidth = 0.25) +
		geom_text(aes(x = stars_loc, label = stars), size = 2.5, color = "#000000") +
		xlim(c(0, enr_max+1)) +
		xlab("") + ylab("") + 
		scale_fill_identity() + scale_color_identity() +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position = "none")

	other  <- ggplot(filt_other, aes(x = enrichment, y = biosample_name)) +
		geom_dots(layout = "swarm", side = "both", shape = 16, point_size = 1.5, linewidth = 0) +
		xlab(x_label) + ylab("") +
		xlim(c(0, enr_max+1)) + 
		scale_fill_identity() + scale_color_identity() +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position = "none")

	# for glufocse: 15
	gr <- cowplot::plot_grid(top, other, ncol = 1, rel_heights = c(1,nrow(filt_top)/5), align = "hv")
	ggsave2(file.path(out_dir, paste0(this_trait, "_enrichment_summary_v2.pdf")), gr, height = 3, width = 3.75)

}

## MAIN
model_id = "scE2G_multiome"
GWAS_benchmark_id <- c("2024_0911_hd_blood", "2024_0829_pancreas")
er_files = file.path("/oak/stanford/groups/engreitz/Users/sheth/GWAS_benchmarking_working/GWAS_E2G_benchmarking/results",
	GWAS_benchmark_id, model_id, "variant_overlap", "enrichmentRecall.thresholded.traitByBiosample.tsv.gz")
out_dir <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0623_CTS_predictions/GWAS_enrichment_specificity"
dir.create(out_dir)
sample_key <- fread("/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0916_global_properties/config/pred_sample_key.tsv")

# header: trait	nVariantsOverlappingEnhancers	biosample	nVariantsTotal	bpEnhancers	nCommonVariantsOverlappingEnhancers	nCommonVariantsTotal
	# recall	enrichment	method	biosampleGroup	traitGroup	SE_log_enr	CI_enr_low	CI_enr_high	p_enr	p_adjust_enr
	# recall_adjust	SE_recall	CI_recall_low	CI_recall_high
er <- lapply(er_files, fread) %>% rbindlist() %>% as_tibble() %>% distinct() %>%
	dplyr::filter(biosampleGroup == FALSE, biosample != "BMMC22_ID2_hi_myeloid_prog") %>%
	mutate(p_adjust_enr = p.adjust(p_enr, method = "bonferroni"),
		significant_enr = ifelse(p_adjust_enr < 0.05, "significant", "not_significant")) %>%
	left_join(sample_key, by = "biosample")
#fwrite(er, file.path(out_dir, "enrichmentRecall.thresholded.combined.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

traits_of_interest <- c("MCH", "MCV", "Lym")
# lapply(traits_of_interest, plot_enrichment_across_biosamples, er, out_dir)

# compare_enrichment_across_biosamples("MCH", 4, er, out_dir)
# compare_enrichment_across_biosamples("MCV", 4, er, out_dir)
#compare_enrichment_across_biosamples("Glucose", 2, er, out_dir) 

plot_trait_summary("MCH", 3, 5, er, out_dir) 
plot_trait_summary("MCV", 3, 5, er, out_dir) 
plot_trait_summary("Glucose", 2, 9, er, out_dir) 