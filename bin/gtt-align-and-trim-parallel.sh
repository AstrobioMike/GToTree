#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
ORANGE='\033[0;33m'
NC='\033[0m'

tmp_dir=$2
faster_alignment=$3
num_muscle_threads=$4
target_gene_suffix=$5

# removing those genomes that need to be removed based on not having enough hits to the target genes
gtt-parse-fasta-by-headers -i ${tmp_dir}/${1}_hits_filtered.tmp -w ${tmp_dir}/sorted_genomes_to_remove.tmp -o ${tmp_dir}/${1}_hits_filtered${target_gene_suffix} --inverse

# aligning
if [ $faster_alignment == 'true' ]; then
    muscle -super5 ${tmp_dir}/${1}_hits_filtered${target_gene_suffix} -output ${tmp_dir}/${1}_aligned.tmp -threads ${num_muscle_threads} &> /dev/null
else
    muscle -align ${tmp_dir}/${1}_hits_filtered${target_gene_suffix} -output ${tmp_dir}/${1}_aligned.tmp -threads ${num_muscle_threads} &> /dev/null
fi

# trimming
trimal -in ${tmp_dir}/${1}_aligned.tmp -out ${tmp_dir}/${1}_trimmed${target_gene_suffix}.tmp -automated1

# removing linewraps:
sed 's/ .*$//' ${tmp_dir}/${1}_trimmed${target_gene_suffix}.tmp | awk '!/^>/ { printf "%s", $0; n="\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' > ${tmp_dir}/${1}_formatted${target_gene_suffix}.tmp

## adding gap-sequences for genomes missing the current gene ##
# finding here which ones have it
grep ">" ${tmp_dir}/${1}_formatted${target_gene_suffix}.tmp | tr -d ">" | sort > ${tmp_dir}/${1}_genomes_with_gene.tmp

# now getting which ones don't have it
comm -23 ${tmp_dir}/final_genomes_from_all_sources.tmp ${tmp_dir}/${1}_genomes_with_gene.tmp | sort > ${tmp_dir}/${1}_needed_gappers.tmp

# creating gap-sequences if needed
if [ -s ${tmp_dir}/${1}_needed_gappers.tmp ]; then

    # making a headers file for when making fasta in a few steps:
    sed 's/^/>/' ${tmp_dir}/${1}_needed_gappers.tmp > ${tmp_dir}/${1}_needed_headers.tmp

    # getting length of the alignment for the current gene:
    aln_length_tmp=$(sed -n '2p' ${tmp_dir}/${1}_formatted${target_gene_suffix}.tmp | wc -c | tr -s " " | cut -f2 -d " ")
    # subtracting 1 for newline characters 
    aln_length_tmp=$(echo "$aln_length_tmp"-1 | bc)
    # making a string of gaps the length of the alignment for those missing it:
    gap_seq=$(printf "%0.s-" $(seq 1 1 $aln_length_tmp))
    # making as many gap sequences as there are genomes missing the current gene:
    num_genomes_to_add=$(wc -l ${tmp_dir}/${1}_needed_gappers.tmp | tr -s " " "\t" | cut -f2)
    for i in $(cat ${tmp_dir}/${1}_needed_gappers.tmp)
    do
        echo "$gap_seq"
    done > ${tmp_dir}/${1}_gaps.tmp

    # making fasta of those genomes missing the current gene:
    paste -d "\n" ${tmp_dir}/${1}_needed_headers.tmp ${tmp_dir}/${1}_gaps.tmp > ${tmp_dir}/${1}_missing_genomes${target_gene_suffix}.tmp
    # catting the genomes missing the current gene together with those that have it
    cat ${tmp_dir}/${1}_formatted${target_gene_suffix}.tmp ${tmp_dir}/${1}_missing_genomes${target_gene_suffix}.tmp > ${tmp_dir}/${1}${target_gene_suffix}.tmp
else
    mv ${tmp_dir}/${1}_formatted${target_gene_suffix}.tmp ${tmp_dir}/${1}${target_gene_suffix}.tmp
fi

## reordering the final fasta of this gene so that all gene sets can be pasted together at end ##
gtt-reorder-fasta -i ${tmp_dir}/${1}${target_gene_suffix}.tmp -w ${tmp_dir}/final_genomes_from_all_sources.tmp -o ${tmp_dir}/${1}_all_aligned${target_gene_suffix}

printf "\n\n\n   --------------------------------------------------------------------------   \n"
printf "\t    Finished aligning and formatting gene-set ${GREEN}$1${NC}.\n"
printf "   --------------------------------------------------------------------------   \n"
