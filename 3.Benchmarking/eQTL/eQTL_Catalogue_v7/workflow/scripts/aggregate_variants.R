suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)})

variant_files = snakemake@input$all_samples %>% as.character() %>% strsplit(" ") %>% unlist()
outfile = snakemake@output$variants

variants_sep = lapply(variant_files, FUN = fread)
variants_all =  rbindlist(variants_sep)
fwrite(variants_all, file=outfile, col.names=TRUE, sep="\t")