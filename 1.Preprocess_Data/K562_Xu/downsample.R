library(Seurat)
library(Signac)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(genomation)
library(tidyverse)
library(Rsamtools)

MyDownsampleCount = function(matrix.count,
                             percentage,
                             seed = 0){
  set.seed(seed)
  matrix.downsample <-
    matrix(rbinom(nrow(matrix.count) * ncol(matrix.count),
                  size = as.matrix(matrix.count),
                  prob = percentage),
           nrow = nrow(matrix.count))
  
  matrix.downsample <- as(matrix.downsample, "sparseMatrix")
  rownames(matrix.downsample) = rownames(matrix.count)
  colnames(matrix.downsample) = colnames(matrix.count)
  
  matrix.downsample
}

MyDownsampleRnaSingle = function(matrix_rna_count,
                                 barcode,
                                 umi_count_per_cell,
                                 output_dir,
                                 seed = 2024){
  matrix_rna_count.filter = matrix_rna_count[,barcode]
  
  if(umi_count_per_cell == "all"){
    percentage = 1
  } else{
    percentage = length(barcode) * umi_count_per_cell / sum(matrix_rna_count.filter)
  }
  
  matrix_rna_count.downsample = 
    MyDownsampleCount(matrix_rna_count.filter,
                      percentage,
                      seed = seed)
  
  matrix_rna_logtp10k.downsample = Seurat::LogNormalize(matrix_rna_count.downsample)
  
  system(paste("mkdir -p",output_dir))
  write.csv(as.data.frame(matrix_rna_logtp10k.downsample),
            gzfile(paste(output_dir,"rna_matrix.csv.gz",sep = "/")),
            quote = F)
}

MyDownsampleRnaBarcode = function(matrix_rna_count,
                                  barcode,
                                  output_dir,
                                  vector_umi_count_per_cell,
                                  seed = 2024){
  
  MyDownsampleRnaSingle(matrix_rna_count = matrix_rna_count,
                        barcode = barcode,
                        umi_count_per_cell = "all",
                        output_dir = paste(output_dir,"all",sep = "/"),
                        seed = seed)
  
  for (umi_count_per_cell.tmp in vector_umi_count_per_cell) {
    
    MyDownsampleRnaSingle(matrix_rna_count = matrix_rna_count,
                          barcode = barcode,
                          umi_count_per_cell = umi_count_per_cell.tmp,
                          output_dir = paste(output_dir,umi_count_per_cell.tmp,sep = "/"),
                          seed = seed)
  }
}

MyDownsampleAtacSingle = function(atac_fragments,
                                  barcode,
                                  fragment_count_per_cell,
                                  output_dir,
                                  seed = 2024){
  atac_fragments.filter = atac_fragments[atac_fragments$V4 %in% barcode,]
  
  if(fragment_count_per_cell == "all"){
    atac_fragments.downsample = atac_fragments.filter
  } else{
    set.seed(seed)
    index.downsample = sample(1:nrow(atac_fragments.filter),length(barcode) * fragment_count_per_cell)
    atac_fragments.downsample = atac_fragments.filter[1:nrow(atac_fragments.filter) %in% index.downsample,]
  }
  
  system(paste("mkdir -p",output_dir))
  write.table(atac_fragments.downsample,
              paste(output_dir,"atac_fragments.tsv",sep = "/"),
              col.names = F,
              row.names = F,
              quote = F,
              sep = "\t")
  
  bgzip(paste(output_dir,"atac_fragments.tsv",sep = "/"),
        paste(output_dir,"atac_fragments.tsv.gz",sep = "/"))
  
  indexTabix(paste(output_dir,"atac_fragments.tsv.gz",sep = "/"),
             format = "bed")
  system(paste("rm -f",paste(output_dir,"atac_fragments.tsv",sep = "/")))
}

MyDownsampleAtacBarcode = function(atac_fragments,
                                   barcode,
                                   output_dir,
                                   vector_fragment_count_per_cell,
                                   seed = 2024){
  
  MyDownsampleAtacSingle(atac_fragments = atac_fragments,
                         barcode = barcode,
                         fragment_count_per_cell = "all",
                         output_dir = paste(output_dir,"all",sep = "/"),
                         seed = seed)
  
  for (fragment_count_per_cell.tmp in vector_fragment_count_per_cell) {
    
    MyDownsampleAtacSingle(atac_fragments = atac_fragments,
                           barcode = barcode,
                           fragment_count_per_cell = fragment_count_per_cell.tmp,
                           output_dir = paste(output_dir,fragment_count_per_cell.tmp,sep = "/"),
                           seed = seed)
  }
}



seurat.data = readRDS("data/processed/1.import/seurat.data.rds")
matrix_rna_count.all = seurat.data@assays$RNA$counts
colnames(matrix_rna_count.all) = seurat.data$barcode


atac_fragments.all = read.delim("/maps/projects/ralab_nnfc-AUDIT/people/lpm537/project/E2G/processed/10x_multiome_230720/ArrayExpress/cellranger_res/EMTAB11264_01/outs/atac_fragments.tsv.gz",
                                comment.char = "#",
                                header = F)
atac_fragments.all = atac_fragments.all[atac_fragments.all$V4 %in% seurat.data$barcode,]


seed = 2024

barcode.all = unname(seurat.data$barcode)
set.seed(seed)
barcode.5000 = sample(barcode.all,5000)
set.seed(seed)
barcode.2000 = sample(barcode.all,2000)
set.seed(seed)
barcode.1000 = sample(barcode.all,1000)
set.seed(seed)
barcode.500 = sample(barcode.all,500)
set.seed(seed)
barcode.200 = sample(barcode.all,200)
set.seed(seed)
barcode.100 = sample(barcode.all,100)




MyDownsampleRnaBarcode(matrix_rna_count = matrix_rna_count.all,
                       barcode = barcode.all,
                       output_dir = "data/processed/5.downsample/all_50_20_10_5_2_100/all/rna/",
                       vector_umi_count_per_cell = c(15000,
                                                     10000,
                                                     5000,
                                                     2000,
                                                     1000,
                                                     500),
                       seed = 2024)
MyDownsampleRnaBarcode(matrix_rna_count = matrix_rna_count.all,
                       barcode = barcode.5000,
                       output_dir = "data/processed/5.downsample/all_50_20_10_5_2_100/5000/rna/",
                       vector_umi_count_per_cell = c(15000,
                                                     10000,
                                                     5000,
                                                     2000,
                                                     1000,
                                                     500),
                       seed = 2024)
MyDownsampleRnaBarcode(matrix_rna_count = matrix_rna_count.all,
                       barcode = barcode.2000,
                       output_dir = "data/processed/5.downsample/all_50_20_10_5_2_100/2000/rna/",
                       vector_umi_count_per_cell = c(15000,
                                                     10000,
                                                     5000,
                                                     2000,
                                                     1000,
                                                     500),
                       seed = 2024)
MyDownsampleRnaBarcode(matrix_rna_count = matrix_rna_count.all,
                       barcode = barcode.1000,
                       output_dir = "data/processed/5.downsample/all_50_20_10_5_2_100/1000/rna/",
                       vector_umi_count_per_cell = c(15000,
                                                     10000,
                                                     5000,
                                                     2000,
                                                     1000,
                                                     500),
                       seed = 2024)
MyDownsampleRnaBarcode(matrix_rna_count = matrix_rna_count.all,
                       barcode = barcode.500,
                       output_dir = "data/processed/5.downsample/all_50_20_10_5_2_100/500/rna/",
                       vector_umi_count_per_cell = c(15000,
                                                     10000,
                                                     5000,
                                                     2000,
                                                     1000,
                                                     500),
                       seed = 2024)
MyDownsampleRnaBarcode(matrix_rna_count = matrix_rna_count.all,
                       barcode = barcode.200,
                       output_dir = "data/processed/5.downsample/all_50_20_10_5_2_100/200/rna/",
                       vector_umi_count_per_cell = c(15000,
                                                     10000,
                                                     5000,
                                                     2000,
                                                     1000,
                                                     500),
                       seed = 2024)
MyDownsampleRnaBarcode(matrix_rna_count = matrix_rna_count.all,
                       barcode = barcode.100,
                       output_dir = "data/processed/5.downsample/all_50_20_10_5_2_100/100/rna/",
                       vector_umi_count_per_cell = c(15000,
                                                     10000,
                                                     5000,
                                                     2000,
                                                     1000,
                                                     500),
                       seed = 2024)










MyDownsampleAtacBarcode(atac_fragments = atac_fragments.all,
                        barcode = barcode.all,
                        output_dir = "data/processed/5.downsample/all_50_20_10_5_2_100/all/atac/",
                        vector_fragment_count_per_cell = c(30000,
                                                           20000,
                                                           10000,
                                                           4000,
                                                           2000,
                                                           1000),
                        seed = 2024)
MyDownsampleAtacBarcode(atac_fragments = atac_fragments.all,
                        barcode = barcode.5000,
                        output_dir = "data/processed/5.downsample/all_50_20_10_5_2_100/5000/atac/",
                        vector_fragment_count_per_cell = c(30000,
                                                           20000,
                                                           10000,
                                                           4000,
                                                           2000,
                                                           1000),
                        seed = 2024)
MyDownsampleAtacBarcode(atac_fragments = atac_fragments.all,
                        barcode = barcode.2000,
                        output_dir = "data/processed/5.downsample/all_50_20_10_5_2_100/2000/atac/",
                        vector_fragment_count_per_cell = c(30000,
                                                           20000,
                                                           10000,
                                                           4000,
                                                           2000,
                                                           1000),
                        seed = 2024)
MyDownsampleAtacBarcode(atac_fragments = atac_fragments.all,
                        barcode = barcode.1000,
                        output_dir = "data/processed/5.downsample/all_50_20_10_5_2_100/1000/atac/",
                        vector_fragment_count_per_cell = c(30000,
                                                           20000,
                                                           10000,
                                                           4000,
                                                           2000,
                                                           1000),
                        seed = 2024)
MyDownsampleAtacBarcode(atac_fragments = atac_fragments.all,
                        barcode = barcode.500,
                        output_dir = "data/processed/5.downsample/all_50_20_10_5_2_100/500/atac/",
                        vector_fragment_count_per_cell = c(30000,
                                                           20000,
                                                           10000,
                                                           4000,
                                                           2000,
                                                           1000),
                        seed = 2024)
MyDownsampleAtacBarcode(atac_fragments = atac_fragments.all,
                        barcode = barcode.200,
                        output_dir = "data/processed/5.downsample/all_50_20_10_5_2_100/200/atac/",
                        vector_fragment_count_per_cell = c(30000,
                                                           20000,
                                                           10000,
                                                           4000,
                                                           2000,
                                                           1000),
                        seed = 2024)
MyDownsampleAtacBarcode(atac_fragments = atac_fragments.all,
                        barcode = barcode.100,
                        output_dir = "data/processed/5.downsample/all_50_20_10_5_2_100/100/atac/",
                        vector_fragment_count_per_cell = c(30000,
                                                           20000,
                                                           10000,
                                                           4000,
                                                           2000,
                                                           1000),
                        seed = 2024)

