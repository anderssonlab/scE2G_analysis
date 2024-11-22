function compute_number_of_unique_predictions() {
    local pred_files=($1)
	local  pred_group=$2
	local temp_out=$3
	mkdir $temp_out	
	current_file=$temp_out/${pred_group}_current_file.tsv

	zcat ${pred_files[1]} > $current_file
	#wc -l $current_file

	n=${#pred_files[@]} 
	for ((idxA=1; idxA<n; idxA++)); do  
		local next_file=${pred_files[idxA]}
		bedtools intersect -a $next_file -b $current_file -f 0.5 -F 0.5 -wa -v >> $current_file

		#wc -l $current_file
	done

	n_predictions=$(cat $current_file | wc -l)

	echo $n_predictions
}

function compute_prediction_overlap() {
	local predA=$1
	local predB=$2

	n_shared=$(bedtools intersect -a $predA -b $predB -f 0.5 -F 0.5 | wc -l)
	n_AnotB=$(bedtools intersect -a $predA -b $predB -f 0.5 -F 0.5 -v | wc -l)
	n_BnotA=$(bedtools intersect -a $predB -b $predA -f 0.5 -F 0.5 -v | wc -l)

	rm $predA
	rm $predB

	echo -e "$n_shared\t$n_AnotB\t$n_BnotA"

}

function make_gene_start_end_files() {
	local sample_key=$1
	local out_dir=$2

	local pred_ids=($(cat $sample_key | sed 1d | cut -f2)) # get all biosamples
	local max=${#pred_ids[@]}   # length of list
	local pred_files=($(cat $sample_key | sed 1d | cut -f3)) # get pred files

	mkdir $out_dir

	for ((idxA=0; idxA<max; idxA++)); do
		echo "Processing sample " $idxA "/" $max ": " ${pred_ids[$idxA]}  
		zcat ${pred_files[$idxA]} | csvtk cut -t -f TargetGene,start,end | sed 1d | sort -k 1,1 -k2,2n | uniq | gzip > $out_dir/${pred_ids[$idxA]}_gene_start_end_threshold.tsv.gz
	done
}

### MAIN
sample_key=$OAK/Users/sheth/scE2G_analysis/2024_0916_global_properties/workflow/n_enhancer_analysis/BMMC_split_scramble_sample_key.tsv
out_dir=$OAK/Users/sheth/scE2G_analysis/2024_0916_global_properties/n_enhancers_BMMC
scratch_dir=/scratch/users/shethm/scE2G_analysis/n_enhancers_BMMC
scripts_dir=$OAK/Users/sheth/scE2G_analysis/2024_0916_global_properties/workflow/n_enhancer_analysis
seed=17
delta=3
n_reps=100

out_file=$out_dir/n_enhancers_${n_reps}r_${seed}s_${delta}d.tsv

#mkdir $out_dir
#mkdir -p $scratch_dir

#make_gene_start_end_files $sample_key $scratch_dir/reformatted_predictions

supergroups=("BMMC5_B" "BMMC5_Dendritic" "BMMC5_Erythroid" "BMMC5_Myeloid" "BMMC5_T")

echo -e "supergroup\trep\tsplit_samples\tscramble_samples\tn_total_split\tn_total_scramble\tn_shared\tn_uniq_split\tn_uniq_scramble" > $out_file

n_sg=${#supergroups[@]} 
for ((idxA=0; idxA<n_sg; idxA++)); do  
	echo "Working on: " ${supergroups[$idxA]}

	for ((rep=0; rep<n_reps; rep++)); do
		R_split=$(Rscript $scripts_dir/get_sample_pred_files.R --supergroup ${supergroups[$idxA]} --dataset "BMMC_split" --seed $seed --delta $delta)
		R_scramble=$(Rscript $scripts_dir/get_sample_pred_files.R --supergroup ${supergroups[$idxA]} --dataset "BMMC_scramble" --seed $seed --delta $delta)

		files_split=$(echo "$R_split" | sed -n '1p')
		files_scramble=$(echo "$R_scramble" | sed -n '1p')

		n_split=$(compute_number_of_unique_predictions "$files_split" "BMMC_split" $scratch_dir | tail -n 1)
		n_scramble=$(compute_number_of_unique_predictions "$files_scramble" "BMMC_scramble" $scratch_dir | tail -n 1)
		overlap=$(compute_prediction_overlap "$scratch_dir/BMMC_split_current_file.tsv" "$scratch_dir/BMMC_scramble_current_file.tsv" | tail -n 1)

		samples_split=$(echo "$R_split" | sed -n '2p')
		samples_scramble=$(echo "$R_scramble" | sed -n '2p')

		echo -e "${supergroups[$idxA]}\t$((rep + 1))\t$samples_split\t$samples_scramble\t$n_split\t$n_scramble\t$overlap" >> $out_file

		seed=$((seed + delta))
	done	 
done