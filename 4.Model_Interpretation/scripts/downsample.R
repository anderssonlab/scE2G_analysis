library(Seurat)
library(Signac)
library(ggplot2)
library(genomation)
library(tidyverse)
library(Rsamtools)
library(Matrix)
library(data.table)

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
	if (file.exists(file.path(output_dir,"rna_matrix.csv.gz"))) {
		message("RNA count matrix already created.")
		for_count = read.csv(file.path(output_dir,"rna_matrix.csv.gz"),
                              row.names = 1,
                              check.names = F)
  		sparse_for_count = Matrix(as.matrix(for_count), sparse = TRUE)
		num_umi <- sum(sparse_for_count)
	} else {
		matrix_rna_count.filter = matrix_rna_count[,barcode]
  
		if(umi_count_per_cell == "all" || umi_count_per_cell == 17961){
			percentage = 1
		} else{
			percentage = length(barcode) * umi_count_per_cell / sum(matrix_rna_count.filter)
		}
		
		matrix_rna_count.downsample = 
			MyDownsampleCount(matrix_rna_count.filter,
							percentage,
							seed = seed)
		
		#matrix_rna_logtp10k.downsample = Seurat::LogNormalize(matrix_rna_count.downsample)
		num_umi = sum(matrix_rna_count.downsample)

		system(paste("mkdir -p",output_dir))
		write.csv(as.data.frame(matrix_rna_count.downsample),
					gzfile(file.path(output_dir,"rna_matrix.csv.gz")),
					quote = F)
	}

	return(num_umi)
  
}

MyDownsampleAtacSingle = function(atac_fragments,
                                  barcode,
                                  fragment_count_per_cell,
                                  output_dir,
                                  seed = 2024){
  
    if (file.exists(file.path(output_dir,"atac_fragments.tsv.gz.tbi")) &
		file.exists(file.path(output_dir,"atac_fragments.tsv.gz")) &
		!file.exists(file.path(output_dir,"atac_fragments.tsv"))) {
			message("ATAC fragment file already created.")
		for_count <-  read.delim(file.path(output_dir,"atac_fragments.tsv.gz"), comment.char = "#", header = F)
		num_frag <- nrow(for_count)
	} else {
		atac_fragments.filter = atac_fragments[atac_fragments$V4 %in% barcode,]

		if(fragment_count_per_cell == "all" || fragment_count_per_cell == 33387){
			atac_fragments.downsample = atac_fragments.filter
		} else{
			set.seed(seed)
			index.downsample = sample(1:nrow(atac_fragments.filter),length(barcode) * fragment_count_per_cell)
			atac_fragments.downsample = atac_fragments.filter[1:nrow(atac_fragments.filter) %in% index.downsample,]
		}
		
		system(paste("mkdir -p",output_dir))
		write.table(atac_fragments.downsample,
					file.path(output_dir,"atac_fragments.tsv"),
					col.names = F,
					row.names = F,
					quote = F,
					sep = "\t")
		num_frag <- nrow(atac_fragments.downsample)

		bgzip(file.path(output_dir,"atac_fragments.tsv"),
				file.path(output_dir,"atac_fragments.tsv.gz"))
		
		indexTabix(file.path(output_dir,"atac_fragments.tsv.gz"),
					format = "bed")
		system(paste("rm -f",file.path(output_dir,"atac_fragments.tsv")))
	}
	return(num_frag)
}


seurat_input = "/oak/stanford/groups/engreitz/Projects/scE2G/data/Xu_et_al_downsampling/input/seurat.data.rds"
atac_input = "/oak/stanford/groups/engreitz/Projects/scE2G/data/Xu_et_al_downsampling/input/atac_fragments.tsv.gz"
out_dir = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2025_0124_new_downsample/scE2G_input"
atac_to_umi_ratio = 10
stats_file <- file.path(out_dir, paste0("stats_ratio_", atac_to_umi_ratio, ".tsv"))

system(paste("mkdir -p",out_dir))

seurat.data = readRDS(seurat_input)
matrix_rna_count.all = seurat.data@assays$RNA$counts
colnames(matrix_rna_count.all) = seurat.data$barcode

atac_fragments.all = read.delim(atac_input,
                                comment.char = "#",
                                header = F)
atac_fragments.all = atac_fragments.all[atac_fragments.all$V4 %in% seurat.data$barcode,]


num_cells = c(100, 200, 500, 1000, 2000, 5000, 7821)
num_frag = c(1000, 2000, 4000, 10000, 20000, 30000, 33387)
max_umi = 17961

num_umi = round(pmin(num_frag / atac_to_umi_ratio, max_umi))

seed_init = 17 * atac_to_umi_ratio

barcodes_all = unname(seurat.data$barcode)

seed = seed_init
counter = 1
stats_list <- vector("list", length(num_cells) * length(num_frag))

for (i in 1:length(num_cells)){ 
	for (j in 1:length(num_frag)) {
		set.seed(seed)

		sample_id = paste0("atac_", num_frag[j], "_rna_", num_umi[j], "_cells_", num_cells[i])
		message(sample_id, " (", counter, "/", length(num_cells) * length(num_frag), ")")
		out_dir_this <- file.path(out_dir, sample_id)

		barcodes_this <- sample(barcodes_all, num_cells[i])

		# RNA matrix
		res_umi <- MyDownsampleRnaSingle(matrix_rna_count = matrix_rna_count.all,
                                 barcode =  barcodes_this,
                                 umi_count_per_cell = num_umi[j],
                                 output_dir = out_dir_this,
                                 seed = seed)

		# ATAC fragments
		res_frag <- MyDownsampleAtacSingle(atac_fragments = atac_fragments.all,
                           barcode = barcodes_this,
                           fragment_count_per_cell = num_frag[j],
                           output_dir = out_dir_this,
                           seed = seed)

		res_row = data.frame(as.list(c(ratio = atac_to_umi_ratio, sample_id = sample_id, num_umi = res_umi, num_frag = res_frag,
			frag_per_cell = num_frag[j], umi_per_cell = num_umi[j], num_cells = num_cells[i])))
		print(res_row)

		stats_list[[counter]] = res_row

		seed <- seed + 3
		counter <- counter + 1
	}
}

stats <- rbindlist(stats_list) %>% as.data.frame()
fwrite(stats, stats_file, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
	
