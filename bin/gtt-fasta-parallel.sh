#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

tmp_dir=$2
hmm_file=$3
num_cpus=$4
hmm_target_genes_total=$5
output_dir=$6

# setting assembly name as filename with no extension
assembly="$(basename ${1%.*})"

# adding assembly to ongoing genomes list
echo $assembly >> ${tmp_dir}/fasta_genomes_list.tmp

num=$((num+1)) # to track progress

## running prodigal to get coding sequences
prodigal -c -q -i $1 -a ${tmp_dir}/${assembly}_genes1.tmp > /dev/null
tr -d '*' < ${tmp_dir}/${assembly}_genes1.tmp > ${tmp_dir}/${assembly}_genes2.tmp

## renaming seqs to have assembly name
gtt-rename-fasta-headers -i ${tmp_dir}/${assembly}_genes2.tmp -w $assembly -o ${tmp_dir}/${assembly}_genes.tmp
  
### running hmm search ###
hmmsearch --cut_ga --cpu $num_cpus --tblout ${tmp_dir}/${assembly}_curr_hmm_hits.tmp $hmm_file ${tmp_dir}/${assembly}_genes.tmp > /dev/null

### calculating % completion and redundancy ###
for SCG in $(cat ${tmp_dir}/uniq_hmm_names.tmp)
do
    grep -w -c "$SCG" ${tmp_dir}/${assembly}_curr_hmm_hits.tmp
done > ${tmp_dir}/${assembly}_uniq_counts.tmp

## adding SCG-hit counts to table
paste <(printf $assembly) <(printf %s "$(cat ${tmp_dir}/${assembly}_uniq_counts.tmp | tr "\n" "\t")") >> ${output_dir}/All_genomes_SCG_hit_counts.tsv

num_SCG_hits=$(awk ' $1 > 0 ' ${tmp_dir}/${assembly}_uniq_counts.tmp | wc -l | tr -s " " | cut -f2 -d " ")

num_SCG_redund=$(awk '{ if ($1 == 0) { print $1 } else { print $1 - 1 } }' ${tmp_dir}/${assembly}_uniq_counts.tmp | awk '{ sum += $1 } END { print sum }')

perc_comp=$(echo "$num_SCG_hits / $hmm_target_genes_total * 100" | bc -l)
perc_comp_rnd=$(printf "%.2f\n" $perc_comp)
perc_redund=$(echo "$num_SCG_redund / $hmm_target_genes_total * 100" | bc -l)
perc_redund_rnd=$(printf "%.2f\n" $perc_redund)

# adding NA for taxid so final table can still have the column and lineage for those that do have them
taxid="NA"

## writing summary info to table ##
printf "$assembly\t$1\t$taxid\t$num_SCG_hits\t$perc_comp_rnd\t$perc_redund_rnd\n" >> ${output_dir}/Fasta_genomes_summary_info.tsv

### Pulling out hits for this genome ###
# looping through ribosomal proteins and pulling out each first hit (hmm results tab is sorted by e-value):
esl-sfetch --index ${tmp_dir}/${assembly}_genes.tmp > /dev/null
        
for SCG in $(cat ${tmp_dir}/uniq_hmm_names.tmp)
do
    grep -w "$SCG" ${tmp_dir}/${assembly}_curr_hmm_hits.tmp | awk '!x[$3]++' | cut -f1 -d " " | esl-sfetch -f ${tmp_dir}/${assembly}_genes.tmp - | sed "s/>.*$/>$assembly/" | sed 's/^Usage.*$//' | sed 's/^To see.*$//' | sed '/^$/d' >> ${tmp_dir}/${SCG}_hits.faa
done

rm -rf ${tmp_dir}/${assembly}_*.tmp ${tmp_dir}/${assembly}_genes.tmp.ssi

printf "    ${GREEN}$assembly${NC} finished.\n"
printf "        Found $num_SCG_hits of the targeted $hmm_target_genes_total genes.\n\n"

