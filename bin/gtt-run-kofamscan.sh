#!/usr/bin/env bash

assembly_id=${1}
genes_file=${2}
target_KOs_table_file=${3}
target_KO_hmms_dir=${4}
tmp_dir=${5}
output_dir=${6}
unique_target_KOs=${7}

curr_genome_output_dir="${output_dir}/KO_searching/individual_genome_results/${assembly_id}/"
mkdir -p ${curr_genome_output_dir}

output_results_table_file="${curr_genome_output_dir}/kofamscan-results.tsv"
tmp_ko_working_dir="${tmp_dir}/kofamscan/${assembly_id}/"
tmp_unique_ko_hits="${tmp_ko_working_dir}/unique-KOs.txt"
output_counts_file="${tmp_ko_working_dir}/KO-counts.txt"

# running scan
exec_annotation -p ${target_KO_hmms_dir} -k ${target_KOs_table_file} --cpu 1 -f mapper --no-report-unannotated --tmp-dir ${tmp_ko_working_dir} -o ${output_results_table_file} ${genes_file}

# moving forward only if there were any hits
if [ -s ${output_results_table_file} ]; then

    # getting all unique KOs with hits in this genome
    cut -f 2 ${output_results_table_file} | sort -u > ${tmp_unique_ko_hits}

    # creating individual fasta files with hits for each
    for ko in $(cat ${tmp_unique_ko_hits}); do

        # getting gene IDs with hits to the current KO
        grep -w ${ko} ${output_results_table_file} | cut -f 1 > ${tmp_ko_working_dir}/${ko}-gene-IDs.txt

        # pulling out seqs for this genome
        gtt-parse-fasta-by-headers -i ${genes_file} -w ${tmp_ko_working_dir}/${ko}-gene-IDs.txt -o ${tmp_ko_working_dir}/${ko}.faa

        # removing gene IDs file
        rm ${tmp_ko_working_dir}/${ko}-gene-IDs.txt

    done

else

    printf "No hits detected.\n" > ${output_results_table_file}

fi

# creating count file that can be stuck together at end
rm -rf ${output_counts_file}

for ko in $(cat ${unique_target_KOs}); do

    grep -w -c ${ko} ${output_results_table_file} >> ${output_counts_file}

done
