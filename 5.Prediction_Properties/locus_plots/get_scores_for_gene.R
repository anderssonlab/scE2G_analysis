
# required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(BiocParallel)
  library(GenomicRanges)
  library(optparse)
})

# register parallel backend
register(MulticoreParam(workers = 4))

# read in arguments
option_list <- list(
	make_option(c("--gene"), type="character", default=NA, help="gene symbol to collect predictions from"),
	make_option(c("--sample_key"), type="character", default="/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0909_locus_plots/reference/sample_key.tsv", help="file with biosample,pred_file"),
	make_option(c("--score_threshold"), type="numeric", default = 0.164, help="min score for thresholded scores"),
	make_option(c("--score_column"), type="character", default = "E2G.Score.qnorm", help="name of score column in prediction files"),
	make_option(c("--out_dir"), type="character", default = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0909_locus_plots/", help="output directory"))

opt = parse_args(OptionParser(option_list=option_list))
gene = opt$gene
sample_key_file = opt$sample_key
e2g_threshold = opt$score_threshold
score_column = opt$score_column
out_dir = file.path(opt$out_dir, gene)

# get list of all samples and files
sample_key = fread(sample_key_file)
e2g_files = sample_key$pred_file
names(e2g_files) = sample_key$biosample

# columns to save
pred_columns <- c("chr", "start", "end", "name", "class", "TargetGene", "TargetGeneTSS", "TargetGeneIsExpressed", "TargetGeneEnsembl_ID",
	"isSelfPromoter",  "CellType", "distance", "normalizedATAC_prom", "ABC.Score", "numCandidateEnhGene", "numTSSEnhGene",
	"numNearbyEnhancers", "ubiqExpressed", "RNA_meanLogNorm", "RNA_pseudobulkTPM", "RNA_percentCellsDetected", "Kendall",
	"ARC.E2G.Score", "E2G.Score", "E2G.Score.qnorm", "E2G.Score.qnorm.ignoreTPM")

# function to load predictions for one sample, extract scores for one gene
extract_scores_gene <- function(file, gene, score_col, threshold = 0) {
  # load predictions and extract scores for desired gene
  pred <- fread(file)
  pred <- filter(pred, TargetGene == gene)
  
  # create GRanges object with scores for all enhancers for desired gene
  scores <- pred %>% 
    filter(class != "promoter") %>% 
    select(chr, start, end, all_of(pred_columns)) %>% 
	mutate(score = !!sym(score_col)) %>%
    filter(score >= threshold)
  
  # create GRanges track from scores data frame
  if (nrow(scores) > 0) {
    scores <- makeGRangesFromDataFrame(scores, keep.extra.columns = TRUE)
  } else {
    scores <- GRanges()
  }
  
  rm(pred)
  invisible(gc())
  
  return(scores)
  
}

# get scores for all  samples
pred_scores <- bplapply(e2g_files, FUN = extract_scores_gene, gene = gene, score_col = score_column,
                        threshold = 0)

# create threshold predictions scores
pred_scores_th <- lapply(pred_scores, FUN = function(x) {
  x[x$score >= e2g_threshold]
})

dir.create(out_dir)

# save pred scores to .rds file (named lists (names = biosamples) of GRanges objects)
saveRDS(pred_scores, file = file.path(out_dir, "all_e2g_scores.rds"))
saveRDS(pred_scores_th , file = file.path(out_dir, "thresholded_e2g_scores.rds"))
