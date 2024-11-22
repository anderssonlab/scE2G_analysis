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

# MyDownsampleRnaBarcode = function(matrix_rna_count,
#                                   barcode,
#                                   output_dir,
#                                   vector_umi_count_per_cell,
#                                   seed = 2024){
  
#   MyDownsampleRnaSingle(matrix_rna_count = matrix_rna_count,
#                         barcode = barcode,
#                         umi_count_per_cell = "all",
#                         output_dir = paste(output_dir,"all",sep = "/"),
#                         seed = seed)
  
#   for (umi_count_per_cell.tmp in vector_umi_count_per_cell) {
#     MyDownsampleRnaSingle(matrix_rna_count = matrix_rna_count,
#                           barcode = barcode,
#                           umi_count_per_cell = umi_count_per_cell.tmp,
#                           output_dir = paste(output_dir,umi_count_per_cell.tmp,sep = "/"),
#                           seed = seed)
#   }
# }

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

# MyDownsampleAtacBarcode = function(atac_fragments,
#                                    barcode,
#                                    output_dir,
#                                    vector_fragment_count_per_cell,
#                                    seed = 2024){
  
#   MyDownsampleAtacSingle(atac_fragments = atac_fragments,
#                          barcode = barcode,
#                          fragment_count_per_cell = "all",
#                          output_dir = paste(output_dir,"all",sep = "/"),
#                          seed = seed)
  
#   for (fragment_count_per_cell.tmp in vector_fragment_count_per_cell) {
    
#     MyDownsampleAtacSingle(atac_fragments = atac_fragments,
#                            barcode = barcode,
#                            fragment_count_per_cell = fragment_count_per_cell.tmp,
#                            output_dir = paste(output_dir,fragment_count_per_cell.tmp,sep = "/"),
#                            seed = seed)
#   }
# }

## MAIN
out_dir = "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0916_global_properties/BMMC_split_input"
dir.create(out_dir)

seurat_input <- ("/oak/stanford/groups/engreitz/Projects/scE2G/data/10x_BMMC_neurips2021/obj.seurat.rds")
#atac_input <- "/oak/stanford/groups/engreitz/Projects/scE2G/data/10x_BMMC_neurips2021/fragments.merge.gz"
obj.seurat <- readRDS(seurat_input)

cell_types <- unique(obj.seurat@meta.data$cell_type.rename)
print(cell_types)

stats_list <- vector("list", length(cell_types))

for (n in seq_along(cell_types)) {
	this_cell_type <- cell_types[n]; message(this_cell_type)
	seurat.subset = subset(obj.seurat, subset = cell_type.rename == this_cell_type)
	
	matrix_rna_count.subset = seurat.subset@assays$RNA$counts
	colnames(matrix_rna_count.subset) = seurat.subset$sample_barcode

	num_cells_total <- length(unique(seurat.subset$sample_barcode)) # number of cells in this cell type
	message("Number of cells in cell type: ", num_cells_total) 

	# calculate numbr of subclusters to make
	if (num_cells_total < 1000) {
		n_clusters = 2
	} else {
		n_clusters = floor(num_cells_total/500)
	}

	atac_input_file <- paste0("/oak/stanford/groups/engreitz/Projects/scE2G/data/10x_BMMC_neurips2021/fragments/", this_cell_type, ".atac_fragments.tsv.gz")
	atac_fragments.this <- read.delim(atac_input_file,
									comment.char = "#",
		 							header = F) # %>%
		# separate_wider_delim(V4, delim = "_", names = c("prefix", "barcode")) %>%
		# mutate(V4 = barcode)
	print(head(atac_fragments.this))

	# get barcodes of cells per sample
	# seed_init = 17
	seed_init = 17
	num_cells = rep(floor(num_cells_total/n_clusters), times = n_clusters)
	num_frag = rep("all", times =  n_clusters)
	num_umi = rep("all", times = n_clusters)
	num_reps =  n_clusters
	barcodes_all = unname(seurat.subset$sample_barcode)
	print(head(barcodes_all))
	print(length(barcodes_all))

	message("Number of subclusters: ", n_clusters, "; cells per subcluster: ", num_cells[1])
	barcodes_per_dimension <- list("vector", n_clusters) # list per cell x frag pair
	seed = seed_init

	cluster_stats_list = vector("list", n_clusters)
	for (i in 1:n_clusters){
		sample_id = paste0(this_cell_type, "_subcluster_", i)
		message(sample_id)
		set.seed(seed)

		barcodes_this <- sample(barcodes_all, num_cells[i])
		print(head(barcodes_this))
		print(length(barcodes_this))
		barcodes_all <- setdiff (barcodes_all, barcodes_this) # remove barcodes used in this rep
		#atac_fragments.this = atac_fragments.all[atac_fragments.all$V4 %in% barcodes_this,]

		out_dir_this <- file.path(out_dir, sample_id)
			
			# RNA matrix
			res_umi <- MyDownsampleRnaSingle(matrix_rna_count = matrix_rna_count.subset,
									barcode =  barcodes_this,
									umi_count_per_cell = num_umi[i],
									output_dir = out_dir_this,
									seed = seed)

			# ATAC fragments
			res_frag <- MyDownsampleAtacSingle(atac_fragments = atac_fragments.this,
							barcode = barcodes_this,
							fragment_count_per_cell = num_frag[i],
							output_dir = out_dir_this,
							seed = seed)
			res_row = data.frame(as.list(c(cell_type = this_cell_type, sample_id = sample_id, num_umi = res_umi, num_frag = res_frag, num_cells = num_cells[i])))
			print(res_row)

			cluster_stats_list[[i]] = res_row
			#seed <- seed + 3
		}

	stats_list[[n]] = rbindlist(cluster_stats_list) %>% as.data.frame()

}

stats <- rbindlist(stats_list) %>% as.data.frame()
fwrite(stats, file.path(out_dir, "stats_v1.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)