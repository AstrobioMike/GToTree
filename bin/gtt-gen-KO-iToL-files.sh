#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
ORANGE='\033[0;33m'
NC='\033[0m'

tmp_dir=$1
output_dir=$2

## setting variable holding whether or not any labels were swapped
if grep -q label <(head -n 1 ${output_dir}/Genomes_summary_info.tsv); then 
    labels_swapped='true'
else
    labels_swapped='false'
fi

curr_target_line=0

for target in $(cat ${tmp_dir}/uniq_ko_targets.tmp)
do

    curr_target_line=$(($curr_target_line + 1))

    target_col=$(($curr_target_line + 2))

    awk -F $'\t' -v col="$target_col" ' $col > 0 { print $1 } ' ${output_dir}/KO_search_results/KO-hit-counts.tsv | tail -n +2 > ${tmp_dir}/Genomes_with_hits_to_${target}.tmp

    if [ -s ${tmp_dir}/Genomes_with_hits_to_${target}.tmp ]; then

        if [ $labels_swapped == 'true' ]; then

            for genome in $(cat ${tmp_dir}/Genomes_with_hits_to_${target}.tmp)
            do
                
                grep -m1 "^$genome" ${output_dir}/Genomes_summary_info.tsv | cut -f 2 

            done > ${tmp_dir}/Genome_labels_with_hits_to_${target}.tmp

            paste ${tmp_dir}/Genomes_with_hits_to_${target}.tmp ${tmp_dir}/Genome_labels_with_hits_to_${target}.tmp > ${tmp_dir}/Genomes_with_hits_to_${target}.tsv

        fi
        
        ## if any, removing those not in final tree before making iToL file
        awk -F $'\t' ' $8 == "No" { print $1 } ' ${output_dir}/Genomes_summary_info.tsv | sort > ${tmp_dir}/sorted_genomes_to_leave_out_of_KO_iToL_files.tmp

        if [ -s ${tmp_dir}/sorted_genomes_to_leave_out_of_KO_iToL_files.tmp ]; then
            comm -23 <( sort ${tmp_dir}/Genomes_with_hits_to_${target}.tmp ) ${tmp_dir}/sorted_genomes_to_leave_out_of_KO_iToL_files.tmp > ${tmp_dir}/genomes_retained_for_${target}_iToL.tmp
        else
            cp ${tmp_dir}/Genomes_with_hits_to_${target}.tmp ${tmp_dir}/genomes_retained_for_${target}_iToL.tmp
        fi

        ## making iToL file for each target KO
        if [ $labels_swapped == 'true' ]; then
            for genome in $(cat ${tmp_dir}/genomes_retained_for_${target}_iToL.tmp)
            do
                grep -m 1 -w "$genome" ${tmp_dir}/Genomes_with_hits_to_${target}.tsv
            done | cut -f 2 > ${tmp_dir}/genomes_for_iToL_for_${target}.tmp

        else
            for genome in $(cat ${tmp_dir}/genomes_retained_for_${target}_iToL.tmp)
            do
                grep -m 1 -w "$genome" ${tmp_dir}/Genomes_with_hits_to_${target}.tmp
            done > ${tmp_dir}/genomes_for_iToL_for_${target}.tmp
        fi

        printf "DATASET_STYLE\nSEPARATOR SPACE\nDATASET_LABEL $target\nCOLOR #0000ff\nDATA\n" > ${output_dir}/KO_search_results/iToL_files/${target}-iToL.txt

        cat <(sed 's/$/ branch node #0000ff 3 normal/' ${tmp_dir}/genomes_for_iToL_for_${target}.tmp) >> ${output_dir}/KO_search_results/iToL_files/${target}-iToL.txt
    
    else
        rm ${tmp_dir}/Genomes_with_hits_to_${target}.tmp

    fi

done
