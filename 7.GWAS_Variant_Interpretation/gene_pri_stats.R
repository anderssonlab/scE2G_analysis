suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
})

out_dir <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0823_GWAS_variant_interpretation"
out_file <- file.path(out_dir, "gene_prioritization_stats.txt")
variant_file = "/oak/stanford/groups/engreitz/Users/sheth/GWAS_benchmarking_working/GWAS_E2G_benchmarking/results/2024_0829_pancreas/variants/filteredGWASVariants.merged.sorted.tsv.gz"
ds <- c("BMMC22_PBMC9", "pancreatic_islets")
input_files <- file.path(out_dir, ds, "gene_prioritization_table_enr_filter.tsv")
E2GRank_lim <- 2
# columns: CredibleSet	trait	TargetGene	bestE2GVariantScore	bestE2GCellType	bestE2GVariantID	nCellTypes	listCellTypes	anyTPMOver1	listTPMs
# E2GRank	bestDNCDistanceToTSS	distanceRank	bestE2GVariantPIP	AnyCoding	AnyPromoter	AnySpliceSite	
# BestOverallSNP	BestSNPOverallNearestGene	isDistalNoncoding	isNearestGeneToBestSNP	sumEnhancersPerGene

df <- lapply(input_files, fread) %>% rbindlist() %>% as_tibble() %>%
	dplyr::filter(isDistalNoncoding %in% c(NA, TRUE))

# INPUT
cs_considered <- df %>% dplyr::select(CredibleSet, trait) %>% distinct()
traits <- unique(cs_considered$trait)
n_traits_cons  <- paste0("Number of  traits  considered: ", length(traits))
n_cs_cons <- paste0("Number of nc cs considered: ", length(unique(cs_considered$CredibleSet)))
n_cs_trait_cons <- paste0("Number of nc cs - trait pairs considered: ", nrow(cs_considered))

n_var <- fread(variant_file) %>%
	dplyr::filter(trait %in% traits, CredibleSet %in% cs_considered$CredibleSet) %>%
	dplyr::select(trait, rsid) %>% distinct()
n_var_cons <- paste0("Number of variants considered: ", length(unique(n_var$rsid)))

# AGGREGATE PRED
df_pred <- df %>% dplyr::filter(E2GRank <= E2GRank_lim)
n_var_pri <- paste0("Number of variants prioritized: ", length(unique(df_pred$bestE2GVariantID)))
n_genes_pri <-  paste0("Number of genes prioritized: ", length(unique(df_pred$TargetGene)))

df_n <- df_pred %>% dplyr::select(CredibleSet, trait) %>% distinct()
n_cs_trait_pri <-  paste0("Number of nc cs-trait pairs prioritized: ", nrow(df_n))
n_cs_pri  <- paste0("Number of nc cs  prioritized: ", length(unique(df_n$CredibleSet)))

df_x <- df_pred %>% dplyr::select(TargetGene, bestE2GVariantID) %>% distinct()
n_pred_total <- paste0("Number of variant-gene predictions: ", nrow(df_x))

df_ct <- df_pred %>% dplyr::select(TargetGene, bestE2GVariantID, bestE2GCellType) %>% distinct()
n_pred_ct <- paste0("Number of variant-gene-best cell type predictions: ", nrow(df_ct))

df_cs <- df_pred %>% dplyr::select(TargetGene, CredibleSet) %>% distinct() %>% nrow()
df_cs_ct <- df_pred %>% dplyr::select(TargetGene, CredibleSet, bestE2GCellType) %>% distinct() %>% nrow()
n_pred_cs <- c(paste0("Number of cs-gene predictions: ", df_cs), 
	paste0("Number of cs-gene-best cell type predictions: ", df_cs_ct))

df_dist <- df_pred %>% dplyr::select(TargetGene, CredibleSet, E2GRank, distanceRank) %>% distinct()
df_n_closest <- df_dist %>% dplyr::filter(distanceRank == 1) %>% nrow()
df_n_dist <- df_dist %>% dplyr::filter(distanceRank != 1) %>% nrow()
df_n_bestE2G_dist <- df_dist %>% dplyr::filter(E2GRank == 1, distanceRank != 1) %>% nrow()
n_dist <- c(paste0("Number of cs-gene pairs linked to closest gene: ", df_n_closest),
	paste0("Number of cs-gene pairs NOT linked to closest gene: ", df_n_dist),
	paste0("Number of cs-gene pairs with E2G Rank 1 not linked to closest gene: ", df_n_bestE2G_dist))

# PRED PROPERTIES 

# PER TRAIT PRED
n_var_per_trait <- n_var %>% 
	group_by(trait) %>%
	summarize(n_variants = n_distinct(rsid))


out_text <- c(n_traits_cons, n_cs_cons, n_cs_trait_cons, n_var_cons,
	n_var_pri, n_genes_pri, n_cs_trait_pri, n_cs_pri, n_pred_total, n_pred_ct,
	n_pred_cs, n_dist)
writeLines(out_text, out_file)