suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)})

variant_file = snakemake@input$variants_filt # columns: gene_id, variant, pip
gene_file = snakemake@params$genes
tissue = snakemake@params$tissue
print(tissue)
outfile = snakemake@output$variants_processed #  columns: chr, start, end, varID_hg38, tissue, gene_hgnc,pip

df = fread(variant_file, sep="\t", header=FALSE)
colnames(df) = c("gene_ensembl", "varID_hg38", "pip")
message("Read variant file.")

# add location columns from var id, in format chr1_108008122_A_G
df = separate(df, varID_hg38, c("chr", "start"), extra = "drop", remove=FALSE, convert=TRUE)
df$end = df$start + 1
message("Added location columns.")

# add/convert to hgncs
genes = fread(gene_file, sep="\t", header=TRUE)
genes = dplyr::select(genes, name, Ensembl_ID)
message("Read gene file.")

colnames(genes) = c("gene_hgnc", "gene_ensembl")
df = inner_join(df, genes, by="gene_ensembl", relationship = "many-to-many") 
message("Mapped gene names.")

# add tissue column
df$tissue = tissue
message("Added tissue column.")

# write to output (columns:  chr, start, end, variant_id, tissue, gene_hgnc, pip)
fwrite(df, file=outfile, col.names=TRUE, sep="\t")
