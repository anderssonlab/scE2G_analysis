suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)
  library(ggplot2)})

variants_file = snakemake@input$variants #  columns: chr, start, end, varID_hg38, tissue, gene_hgnc,gene_ensembl, pip
outfile = snakemake@output$metadata
outplot=snakemake@output$plot

variants = fread(variants_file, sep="\t", header=TRUE)

summ = distinct(variants) %>% 
				dplyr::group_by(tissue) %>%
				summarize(n_variants=n()) %>%
				arrange(n_variants)
summ$tissue = factor(summ$tissue, levels=summ$tissue, ordered=TRUE)

# plot
g=ggplot(data=summ, aes(x=tissue, y=n_variants)) +
	geom_bar(stat="identity") +
	geom_text(aes(y=n_variants+150, label=n_variants), size=2.5) +
	labs(x="", y="Number of eQTLs above PIP threshold") +
	theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8)) +
	coord_flip() 

ggsave(outplot, g, width=6, height=8)
fwrite(summ, file=outfile, col.names=TRUE, sep="\t")