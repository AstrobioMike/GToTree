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

assembly="$(basename ${1%.*})"

rm -rf ${output_dir}/Genbank_files_with_no_CDSs.txt # deleting if file exists

# adding assembly to ongoing genomes list
echo $assembly >> ${tmp_dir}/genbank_genomes_list.tmp

# storing more info about the assembly if it's present in the genbank file:
# checking for organism:
if grep -q "ORGANISM" $1; then 
    org_name=$(grep -m1 "ORGANISM" $1 | tr -s " " | cut -f3- -d " " | tr "[ ./\\]" "_" | tr -s "_")
else
    org_name="NA"
fi

if grep -q "taxon" $1; then
    taxid=$(grep -m1 "taxon" $1 | cut -f2 -d ":" | tr -d '"')
else
    taxid="NA"
fi

# extracting AA coding sequences from genbank file
gtt-genbank-to-AA-seqs -i $1 -o ${tmp_dir}/${assembly}_genes2.tmp

# checking that the file had CDS annotations
if [ ! -s ${tmp_dir}/${assembly}_genes2.tmp ]; then

    printf "     ${RED}******************************* ${NC}NOTICE ${RED}*******************************${NC}  \n"
    printf "\t  $assembly doesn't appear to have CDS annotations, so we\n"
    printf "\t  are identifying coding sequences with prodigal.\n\n"

    printf "\t    Reported in \"${output_dir}/Genbank_files_with_no_CDSs.txt\".\n"
    printf "     ${RED}**********************************************************************${NC}  \n\n"

    echo "$1" >> ${output_dir}/Genbank_files_with_no_CDSs.txt
    rm -rf ${tmp_dir}/${assembly}_genes2.tmp

    # pulling out full nucleotide fasta from genbank file
    gtt-genbank-to-fasta -i $1 -o ${tmp_dir}/${assembly}_fasta.tmp

    printf "      Getting coding seqs...\n\n"

    # running prodigal
    prodigal -c -q -i ${tmp_dir}/${assembly}_fasta.tmp -a ${tmp_dir}/${assembly}_genes1.tmp > /dev/null
    tr -d '*' < ${tmp_dir}/${assembly}_genes1.tmp > ${tmp_dir}/${assembly}_genes2.tmp

fi

## renaming seqs to have assembly name
gtt-rename-fasta-headers -i ${tmp_dir}/${assembly}_genes2.tmp -w $assembly -o ${tmp_dir}/${assembly}_genes.tmp

printf "   ${GREEN}$assembly${NC}\n"
printf "      Performing HMM search...\n"
  
### running hmm search ###
hmmsearch --cut_ga --cpu $num_cpus --tblout ${tmp_dir}/${assembly}_curr_hmm_hits.tmp $hmm_file ${tmp_dir}/${assembly}_genes.tmp > /dev/null

### calculating % completion and redundancy ###
for SCG in $(cat ${tmp_dir}/uniq_hmm_names.tmp)
do
    grep -w -c "$SCG" ${tmp_dir}/${assembly}_curr_hmm_hits.tmp
done > ${tmp_dir}/${assembly}_uniq_counts.tmp

num_SCG_hits=$(awk ' $1 > 0 ' ${tmp_dir}/${assembly}_uniq_counts.tmp | wc -l | tr -s " " | cut -f2 -d " ")

printf "        Found $num_SCG_hits of the targeted $hmm_target_genes_total.\n\n"

num_SCG_redund=$(awk '{ if ($1 == 0) { print $1 } else { print $1 - 1 } }' ${tmp_dir}/${assembly}_uniq_counts.tmp | awk '{ sum += $1 } END { print sum }')

perc_comp=$(echo "$num_SCG_hits / $hmm_target_genes_total * 100" | bc -l)
perc_comp_rnd=$(printf "%.2f\n" $perc_comp)
perc_redund=$(echo "$num_SCG_redund / $hmm_target_genes_total * 100" | bc -l)
perc_redund_rnd=$(printf "%.2f\n" $perc_redund)

## writing summary info to table ##
printf "$assembly\t$1\t$taxid\t$org_name\t$num_SCG_hits\t$perc_comp_rnd\t$perc_redund_rnd\n" >> ${output_dir}/Genbank_genomes_summary_info.tsv

### Pulling out hits for this genome ###
# looping through ribosomal proteins and pulling out each first hit (hmm results tab is sorted by e-value):
esl-sfetch --index ${tmp_dir}/${assembly}_genes.tmp > /dev/null
        
for SCG in $(cat ${tmp_dir}/uniq_hmm_names.tmp)
do
    grep -w "$SCG" ${tmp_dir}/${assembly}_curr_hmm_hits.tmp | awk '!x[$3]++' | cut -f1 -d " " | esl-sfetch -f ${tmp_dir}/${assembly}_genes.tmp - | sed "s/>.*$/>$assembly/" | sed 's/^Usage.*$//' | sed 's/^To see.*$//' | sed '/^$/d' >> ${tmp_dir}/${SCG}_hits.faa
done

rm -rf ${tmp_dir}/${assembly}_*.tmp ${tmp_dir}/${assembly}_genes.tmp.ssi
