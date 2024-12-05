sample_key=$OAK/Users/sheth/scE2G_analysis/2024_0916_global_properties/config/pred_sample_key_basic.tsv
out_dir=$OAK/Users/sheth/scE2G_analysis/2024_0916_global_properties/correlation_across_clusters_threshold
scratch_dir=/scratch/users/shethm/scE2G_analysis/pairwise_comparisons_correlation_threshold

# sample_key=$OAK/Users/sheth/scE2G_analysis/2024_0916_global_properties/config/pred_sample_key_BMMC_split.tsv
# out_dir=$OAK/Users/sheth/scE2G_analysis/2024_0916_global_properties/correlation_across_clusters_threshold_BMMC_split
# scratch_dir=/scratch/users/shethm/scE2G_analysis/pairwise_comparisons_correlation_threshold_BMMC_split

scripts_dir=$OAK/Users/sheth/scE2G_analysis/2024_0916_global_properties/workflow/scripts
out_file=$out_dir/correlation_across_clusters.tsv

mkdir $out_dir
mkdir -p $scratch_dir
#echo -e "biosampleA\tbiosampleB\tnSharedPredAwB\tnTotalPredA\tnTotalPredB\tnAnotB\tnBnotA\tSpearman\tPearson\tPearson_log1p" > $out_file # header

pred_ids=($(cat $sample_key | sed 1d | cut -f1)) # get all biosamples
max=${#pred_ids[@]}   # length of list
pred_files=($(cat $sample_key | sed 1d | cut -f2)) # get pred files

# get file with columns TargetGene,start,end for thresholded predictions
mkdir $scratch_dir/reformatted_thresholded_predictions
# for ((idxA=0; idxA<max; idxA++)); do
# 	echo "Processing file: " ${pred_files[$idxA]} 
# 	zcat ${pred_files[$idxA]} | csvtk cut -t -f TargetGene,start,end,E2G.Score.qnorm | sed 1d | sort -k 1,1 -k2,2n | uniq | gzip > $scratch_dir/reformatted_thresholded_predictions/${pred_ids[$idxA]}_chr_start_end_score_thresholded.tsv.gz
# done

# intersect all pairs and correlate
for ((idxA=0; idxA<max; idxA++)); do  
	echo "idxA = " $idxA " out of " $max
	temp_int_dir=$scratch_dir/temporary_intersections_idx_$idxA
	mkdir $temp_int_dir
	for ((idxB=idxA; idxB<max; idxB++)); do

	A=$scratch_dir/reformatted_thresholded_predictions/${pred_ids[$idxA]}_chr_start_end_score_thresholded.tsv.gz
	B=$scratch_dir/reformatted_thresholded_predictions/${pred_ids[$idxB]}_chr_start_end_score_thresholded.tsv.gz
	nA=$(zcat $A | wc -l) # number of predicitons for A
	nB=$(zcat $B | wc -l) # number of predictions for B

	temp_int=$temp_int_dir/${pred_ids[$idxA]}_${pred_ids[$idxB]}_intersect.tsv.gz
	temp_anotb=$temp_int_dir/${pred_ids[$idxA]}_${pred_ids[$idxB]}_anotb.tsv
	temp_bnota=$temp_int_dir/${pred_ids[$idxA]}_${pred_ids[$idxB]}_bnota.tsv

	bedtools intersect -a $A -b $B -f 0.5 -F 0.5 -wa -wb | gzip > $temp_int
	bedtools intersect -a $A -b $B -f 0.5 -F 0.5 -v -wa > $temp_anotb
	bedtools intersect -a $B -b $A -f 0.5 -F 0.5 -v -wa > $temp_bnota

	n_shared=$(zcat $temp_int | wc -l) # number of shared predictions
	AnotB=$(bedtools intersect -a $A -b $B -f 0.5 -F 0.5 -v -wa | uniq | wc -l)
	BnotA=$(bedtools intersect -a $B -b $A -f 0.5 -F 0.5 -v -wa | uniq | wc -l)

	corr_out=$(Rscript $scripts_dir/compute_e2g_score_correlation.R --intersection_file $temp_int --a_not_b $temp_anotb --b_not_a $temp_bnota)
	spearman=$(echo $corr_out | awk  '{print $1}')
	pearson=$(echo $corr_out | awk  '{print $2}')
	pearson_log1p=$(echo $corr_out | awk  '{print $3}')
	
	echo -e "${pred_ids[$idxA]}\t${pred_ids[$idxB]}\t$n_shared\t$nA\t$nB\t$AnotB\t$BnotA\t$spearman\t$pearson\t$pearson_log1p" >> $out_file
	echo -e "${pred_ids[$idxA]}\t${pred_ids[$idxB]}\t$n_shared\t$nA\t$nB\t$AnotB\t$BnotA\t$spearman\t$pearson\t$pearson_log1p" 
  done

  rm -rf $temp_int_dir

done
