#!/usr/bin/env bash

all_assembly_ids=${1}
tmp_dir=${2}
output_dir=${3}
unique_target_KOs=${4}

KO_hits_fasta_output_dir="${output_dir}/KO_searching/KO_hit_seqs/"
mkdir -p ${KO_hits_fasta_output_dir}

# combining all fasta files for each individual KO
for ko in $(cat ${unique_target_KOs}); do

    find ${tmp_dir}/kofamscan/ -name ${ko}.faa -exec cat {} \; > ${KO_hits_fasta_output_dir}/${ko}.faa

    # removing if there were none
    if [ ! -s ${KO_hits_fasta_output_dir}/${ko}.faa ]; then

        rm ${KO_hits_fasta_output_dir}/${ko}.faa

    fi

done

# combining counts into one table
building_counts_tab="${tmp_dir}/kofamscan/KO-counts-building.tmp"
final_counts_tab="${output_dir}/KO_searching/KO-counts.tsv"


# starting first row
cat <( printf "KO_ID\n" ) ${unique_target_KOs} > ${building_counts_tab}

# looping through assemblies and adding them
for assembly_id in $(cat ${all_assembly_ids}); do

    paste ${building_counts_tab} <( cat <( printf "${assembly_id}\n" ) ${tmp_dir}/kofamscan/${assembly_id}/KO-counts.txt ) > ${building_counts_tab}.2
    mv ${building_counts_tab}.2 ${building_counts_tab}

done

mv ${building_counts_tab} ${final_counts_tab}
