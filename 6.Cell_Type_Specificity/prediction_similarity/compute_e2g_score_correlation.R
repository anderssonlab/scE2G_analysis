suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(ggdist)
  library(cowplot)
  library(optparse)
})

option_list <- list(
    make_option(c("--intersection_file"), type="character", default=NA, help="input file"),
    make_option(c("--a_not_b"), type="character"),
    make_option(c("--b_not_a"), type="character"))
  
opt = parse_args(OptionParser(option_list=option_list))

int <- fread(opt$intersection_file, header = FALSE) %>%
	setNames(c("geneA", "startA", "endA", "scoreA", "geneB", "startB", "endB", "scoreB")) %>%
	dplyr::select(scoreA, scoreB)

if (file.info(opt$a_not_b)$size > 0) {
	a_not_b <- fread(opt$a_not_b, header= FALSE)%>%
		setNames(c("geneA", "startA", "endA", "scoreA")) %>%
		dplyr::select(scoreA) %>%
		mutate(scoreB = 0)
} else {
	a_not_b <- data.frame(scoreA = numeric(0), scoreB = numeric(0))
}

if (file.info(opt$b_not_a)$size > 0) {
	b_not_a <- fread(opt$b_not_a, header= FALSE) %>%
		setNames(c("geneB", "startB", "endB", "scoreB")) %>%
		dplyr::select(scoreB) %>%
		mutate(scoreA = 0)
} else {
	b_not_a <- data.frame(scoreA = numeric(0), scoreB = numeric(0))
}

int_all <- rbind(int, a_not_b, b_not_a)

spearman <- cor(int_all$scoreA, int_all$scoreB, method = "spearman")
pearson <- cor(int_all$scoreA, int_all$scoreB, method = "pearson")
log1p_pearson <- cor(log1p(int_all$scoreA), log1p(int_all$scoreB), method = "pearson")

cat(spearman, pearson, log1p_pearson)
