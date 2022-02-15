#!/usr/bin/env bash

assembly_id=${1}
genes_file=${2}
gene_count=${3}
num_cpus=${4}
tmp_dir=${5}
output_dir=${6}

hmmsearch_output_file="${output_dir}/Pfam_search_results/individual_genome_results/${assembly_id}_hmmsearch.txt"

hmmsearch --cut_ga --cpu $num_cpus --tblout ${hmmsearch_output_file} ${tmp_dir}/all_pfam_targets.hmm ${genes_file} > /dev/null

### getting counts of each target in this genome
for target in $(cat ${tmp_dir}/actual_pfam_targets.tmp)
do
    grep -w ${target} ${hmmsearch_output_file} | wc -l | sed 's/^ *//' >> ${tmp_dir}/${assembly_id}_hit_counts.tmp
done

### writing results to main output file
paste <( printf "${assembly_id}\t${gene_count}" ) <(printf %s "$(cat ${tmp_dir}/${assembly_id}_hit_counts.tmp | tr "\n" "\t" | sed 's/\t$/\n/')" ) >> ${output_dir}/Pfam_search_results/Pfam-hit-counts.tsv

### Pulling out hits to additional pfam targets for this genome ###
for target in $(cat ${tmp_dir}/actual_pfam_targets.tmp)
do
    if grep -w -q "$target" ${hmmsearch_output_file}; then

        grep -w "$target" ${hmmsearch_output_file} | cut -f 1 -d " " >> ${tmp_dir}/${assembly_id}_${target}_genes_of_int.tmp

        for gene in $(cat ${tmp_dir}/${assembly_id}_${target}_genes_of_int.tmp)
        do
            echo $gene | esl-sfetch -f ${tmp_dir}/${assembly_id}_genes.tmp -
        done >> ${tmp_dir}/${assembly_id}_${target}_genes1.tmp

        gtt-append-fasta-headers -i ${tmp_dir}/${assembly_id}_${target}_genes1.tmp -w ${assembly_id}_${target} -o ${tmp_dir}/${assembly_id}_${target}_genes.tmp
    
        # adding to fasta of that target holding all genomes
        cat ${tmp_dir}/${assembly_id}_${target}_genes.tmp >> ${output_dir}/Pfam_search_results/Pfam_hit_seqs/${target}-hits.faa
    fi

done
