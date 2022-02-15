#!/usr/bin/env bash

target_KOs_file=$1
output_dir=$2

full_KO_list_file="${KO_data_dir}/ko_list"
full_KO_hmms_dir="${KO_data_dir}/profiles"

sub_KO_list_file="${output_dir}/KO_search_results/target-KOs.tsv"
sub_KO_hmms_dir="${output_dir}/KO_search_results/target_KO_profiles/"

# making target ko_list file and copying target hmms to working area
head -n 1 ${full_KO_list_file} > ${sub_KO_list_file}

for ko in $(cat ${target_KOs_file}); do
    grep -m 1 -w "^${ko}" ${full_KO_list_file} >> ${sub_KO_list_file}
    cp ${full_KO_hmms_dir}/${ko}.hmm ${sub_KO_hmms_dir}
done
