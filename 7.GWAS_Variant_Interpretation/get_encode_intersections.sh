# bedtools intersect a variant wiht all 14k encode-re2g predictions and list resulting predictions ordered by score

input_id=RASD1_variant
input_bed=$OAK/Users/sheth/scE2G_analysis/2024_0823_GWAS_variant_interpretation/config/RASD1_var_input.bed
e2g_files=$OAK/Users/sheth/scE2G_analysis/2024_0823_GWAS_variant_interpretation/config/encode_re2g_thresholded_pred_paths.tsv
out_dir=$OAK/Users/sheth/scE2G_analysis/2024_0823_GWAS_variant_interpretation/encode_re2g_intersections
random_pred_for_header=/oak/stanford/groups/engreitz/Projects/Benchmarking/Revisions/Predictors/ENCODE-rE2G/dhs_only/thresholded_predictions/encode_e2g_predictions_neuronal_stem_cell_ENCSR278FVO_thresholded_predictions.tsv.gz

#mkdir $out_dir

out_file=$out_dir/$input_id.tsv
pred_header=$(zcat $random_pred_for_header | head -1 |  tr -d "#")
var_header="chr_var\tstart_var\tend_var\tvariant_id\tcs_loc_hg19" # start,end = position,position+1
echo -e $pred_header $var_header | tr ' ' \\t > $out_file

i=1
n_lines=$(cat $e2g_files | wc -l)

while IFS=$'\t' read -r pred_file
do
	if [ `expr $i % 100` -eq 0 ]
	then
		echo "Intersecting  file $i out of $n_lines"
	fi

	bedtools intersect -a $pred_file -b $input_bed -wa -wb >> $out_file
    (( i += 1))

done < "$e2g_files"