library(Seurat)
library(Signac)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(genomation)
library(tidyverse)
library(Rsamtools)
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
  matrix_rna_count.filter = matrix_rna_count[,barcode]
  
  if(umi_count_per_cell == "all"){
    percentage = 1
  } else{
    percentage = length(barcode) * umi_count_per_cell / sum(matrix_rna_count.filter)
  }
	num_umi = sum(matrix_rna_count.filter)
  
  matrix_rna_count.downsample = 
    MyDownsampleCount(matrix_rna_count.filter,
                      percentage,
                      seed = seed)
  
  matrix_rna_logtp10k.downsample = Seurat::LogNormalize(matrix_rna_count.downsample)
  
  system(paste("mkdir -p",output_dir))
  write.csv(as.data.frame(matrix_rna_logtp10k.downsample),
            gzfile(paste(output_dir,"rna_matrix.csv.gz",sep = "/")),
            quote = F)

	return(num_umi)
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
  num_frag <- nrow(atac_fragments.downsample)

  bgzip(paste(output_dir,"atac_fragments.tsv",sep = "/"),
        paste(output_dir,"atac_fragments.tsv.gz",sep = "/"))
  
  indexTabix(paste(output_dir,"atac_fragments.tsv.gz",sep = "/"),
             format = "bed")
  system(paste("rm -f",paste(output_dir,"atac_fragments.tsv",sep = "/")))

  return(num_frag)
}

## MAIN
all_supergroups = c("B", "T", "Myeloid", "Dendritic", "Erythroid")
#all_n_clusters = c(3, 4, 4, 2, 3)
out_dir = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0916_global_properties/BMMC_scramble_500_input"
dir.create(out_dir)

stats_list <- vector("list", length(all_supergroups))
for (n in seq_along(all_supergroups)) {
	supergroup = all_supergroups[n]; message(supergroup)
	#n_clusters = all_n_clusters[n]

	seurat_input = paste0("/oak/stanford/groups/engreitz/Projects/scE2G/data/10x_BMMC_neurips2021_5_supergroups/seurat_obj/obj.seurat.", supergroup, ".rds")
	atac_input = paste0("/oak/stanford/groups/engreitz/Projects/scE2G/data/10x_BMMC_neurips2021_5_supergroups/fragments/", supergroup, ".atac_fragments.tsv.gz")

	seurat.data = readRDS(seurat_input)
	matrix_rna_count.all = seurat.data@assays$RNA$counts
	colnames(matrix_rna_count.all) = seurat.data$sample_barcode

	atac_fragments.all = read.delim(atac_input,
									comment.char = "#",
									header = F)
	atac_fragments.all = atac_fragments.all[atac_fragments.all$V4 %in% seurat.data$sample_barcode,] # redundant

	num_cells_total <- length(unique(seurat.data$barcode))
	message("Number of cells in supergroup: ", num_cells_total) # should be 9752 for B

	# calculate numbr of subclusters to make
	if (num_cells_total < 1000) {
		n_clusters = 2
	} else {
		n_clusters = floor(num_cells_total/500)
	}

	# get barcodes of cells per sample
	seed_init = 19
	num_cells = rep(floor(num_cells_total/n_clusters), times = n_clusters)
	num_frag = rep("all", times =  n_clusters)
	num_umi = rep("all", times = n_clusters)
	num_reps =  n_clusters
	barcodes_all = unname(seurat.data$sample_barcode)
	print(length(barcodes_all))

	message("Number of subclusters: ", n_clusters, "; cells per subcluster: ", num_cells[1])
	barcodes_per_dimension <- list("vector", n_clusters) # list per cell x frag pair
	seed = seed_init

	cluster_stats_list = vector("list", n_clusters)
	for (i in 1:n_clusters){
		sample_id = paste0("supergroup_", supergroup, "_subcluster_", i)
		message(sample_id)
		set.seed(seed)

		barcodes_this <- sample(barcodes_all, num_cells[i])
		barcodes_all <- setdiff (barcodes_all, barcodes_this) # remove barcodes used in this rep
		#atac_fragments.this = atac_fragments.all[atac_fragments.all$V4 %in% barcodes_this,]

		out_dir_this <- file.path(out_dir, sample_id)
			
			# RNA matrix
			res_umi <- MyDownsampleRnaSingle(matrix_rna_count = matrix_rna_count.all,
									barcode =  barcodes_this,
									umi_count_per_cell = num_umi[i],
									output_dir = out_dir_this,
									seed = seed)

			# ATAC fragments
			res_frag <- MyDownsampleAtacSingle(atac_fragments = atac_fragments.all,
							barcode = barcodes_this,
							fragment_count_per_cell = num_frag[i],
							output_dir = out_dir_this,
							seed = seed)
			res_row = data.frame(as.list(c(cell_type = supergroup, sample_id = sample_id, num_umi = res_umi, num_frag = res_frag, num_cells = num_cells[i])))
			print(res_row)

			cluster_stats_list[[i]] = res_row
			#seed <- seed + 3
		}

	stats_list[[n]] = rbindlist(cluster_stats_list) %>% as.data.frame()

}

stats <- rbindlist(stats_list) %>% as.data.frame()
fwrite(stats, file.path(out_dir, "stats_v1.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)