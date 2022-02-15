#!/usr/bin/env bash

all_assembly_ids=${1}
tmp_dir=${2}
output_dir=${3}
unique_target_KOs=${4}

KO_output_dir="${output_dir}/KO_search_results"
KO_hits_fasta_output_dir="${KO_output_dir}/KO_hit_seqs/"

mkdir -p ${KO_hits_fasta_output_dir}

# combining all fasta files for each individual KO
for ko in $(cat ${unique_target_KOs}); do

    find ${tmp_dir}/kofamscan/ -name ${ko}.faa -exec cat {} \; > ${KO_hits_fasta_output_dir}/${ko}-hits.faa

    # removing if there were none
    if [ ! -s ${KO_hits_fasta_output_dir}/${ko}-hits.faa ]; then

        rm ${KO_hits_fasta_output_dir}/${ko}-hits.faa

    fi

done

# combining counts into one table
final_counts_tab="${KO_output_dir}/KO-hit-counts.tsv"


# starting first row
# cat <( printf "KO_ID\n" ) ${unique_target_KOs} > ${building_counts_tab}
paste <( printf "assembly_id\ttotal_gene_count" ) <( tr "\n" "\t" < ${unique_target_KOs} | sed 's/\t$/\n/' ) > ${final_counts_tab}

# looping through assemblies and adding them
for assembly_id in $(cat ${all_assembly_ids}); do

    cat ${tmp_dir}/kofamscan/${assembly_id}/KO-counts.txt >> ${final_counts_tab}

done
