#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

tmp_dir=$(cat temp_dir_name.tmp)
hmm_file=$(cat hmm_file_path.tmp)
fasta_genomes_total=$(cat fasta_genomes_total.tmp)
num_cpus=$(cat num_cpus.tmp)
hmm_target_genes_total=$(cat hmm_target_genes_total.tmp)

# setting assembly name as filename with no extension
assembly="$(basename ${1%.*})"

# adding assembly to ongoing genomes list
echo $assembly >> ${tmp_dir}/fasta_genomes_list.tmp

num=$((num+1)) # to track progress

printf "  $assembly\n"
printf "      Getting coding seqs...\n\n"

## running prodigal to get coding sequences
prodigal -c -q -i $1 -a ${assembly}_genes1.tmp > /dev/null
tr -d '*' < ${assembly}_genes1.tmp > ${assembly}_genes2.tmp

## renaming seqs to have assembly name
gtt-rename-fasta-headers -i ${assembly}_genes2.tmp -w $assembly -o ${assembly}_genes.tmp

printf "      Performing HMM search...\n"
  
### running hmm search ###
hmmsearch --cut_ga --cpu $num_cpus --tblout ${assembly}_curr_hmm_hits.tmp $hmm_file ${assembly}_genes.tmp > /dev/null

### calculating % completion and redundancy ###
for SCG in $(cat ${tmp_dir}/uniq_hmm_names.tmp)
do
    grep -w -c "$SCG" ${assembly}_curr_hmm_hits.tmp
done > ${assembly}_uniq_counts.tmp

num_SCG_hits=$(awk ' $1 > 0 ' ${assembly}_uniq_counts.tmp | wc -l | tr -s " " | cut -f2 -d " ")

printf "        Found $num_SCG_hits of the targeted $hmm_target_genes_total.\n\n"

num_SCG_redund=$(awk '{ if ($1 == 0) { print $1 } else { print $1 - 1 } }' ${assembly}_uniq_counts.tmp | awk '{ sum += $1 } END { print sum }')

perc_comp=$(echo "$num_SCG_hits / $hmm_target_genes_total * 100" | bc -l)
perc_comp_rnd=$(printf "%.2f\n" $perc_comp)
perc_redund=$(echo "$num_SCG_redund / $hmm_target_genes_total * 100" | bc -l)
perc_redund_rnd=$(printf "%.2f\n" $perc_redund)

# adding NA for taxid so final table can still have the column and lineage for those that do have them
taxid="NA"

## writing summary info to table ##
printf "$assembly\t$1\t$taxid\t$num_SCG_hits\t$perc_comp_rnd\t$perc_redund_rnd\n" >> Fasta_genomes_summary_info.tsv

### Pulling out hits for this genome ###
# making fasta file searchable to pull out the hits (Easel 0.45h June 2018)
esl-sfetch --index ${assembly}_genes.tmp > /dev/null

# looping through ribosomal proteins and pulling out each first hit (hmm results tab is sorted by e-value):
# done as a separate loop just for clarity (well, in hopes of clarity)
for SCG in $(cat ${tmp_dir}/uniq_hmm_names.tmp)

do
  grep -w "$SCG" ${assembly}_curr_hmm_hits.tmp | awk '!x[$3]++' | cut -f1 -d " " | esl-sfetch -f ${assembly}_genes.tmp - | sed "s/>.*$/>$assembly/" | sed 's/^Usage.*$//' | sed 's/^To see.*$//' | sed '/^$/d' >> ${tmp_dir}/${SCG}_hits.faa
done

rm ${assembly}_*.tmp
rm ${assembly}_*.tmp.ssi

