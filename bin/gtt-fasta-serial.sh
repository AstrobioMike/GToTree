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

# looping through the lines of the provided [-f] file (this loop operates on one genome at a time)
while IFS=$'\t' read -r -a file
do

    # setting assembly name as filename with no extension
    assembly="$(basename ${file%.*})"

    # adding assembly to ongoing genomes list
    echo $assembly >> ${tmp_dir}/fasta_genomes_list.tmp

    num=$((num+1)) # to track progress

    printf "   --------------------------------------------------------------------------   \n"
    printf "\tOn assembly ${GREEN}$assembly${NC}; Number $num of $fasta_genomes_total total.\n"
    printf "   --------------------------------------------------------------------------   \n\n"

    printf "      Getting coding seqs...\n\n"

    ## running prodigal to get coding sequences
    prodigal -c -q -i $file -a ${assembly}_genes1.tmp > /dev/null
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
    printf "$assembly\t$file\t$taxid\t$num_SCG_hits\t$perc_comp_rnd\t$perc_redund_rnd\n" >> Fasta_genomes_summary_info.tsv


    ### Pulling out hits for this genome ###
        # this was faster with esl-sfetch, but can't figure out how to install that with conda and i don't think it's too bad without it
        # but when i want to improve efficiency, this is a good place to start, it's a tad excessive at the moment
    # looping through ribosomal proteins and pulling out each first hit (hmm results tab is sorted by e-value):
        
    for SCG in $(cat ${tmp_dir}/uniq_hmm_names.tmp)
    do
        grep -w -m1 "$SCG" ${assembly}_curr_hmm_hits.tmp | awk '!x[$3]++' | cut -f1 -d " " > ${assembly}_${SCG}_curr_wanted_id.tmp
        gtt-parse-fasta-by-headers -i ${assembly}_genes.tmp -w ${assembly}_${SCG}_curr_wanted_id.tmp -o ${assembly}_${SCG}_hit.tmp
        sed 's/\(.*\)_.*/\1/' ${assembly}_${SCG}_hit.tmp >> ${tmp_dir}/${SCG}_hits.faa
    done

    rm -rf ${assembly}_*.tmp

done < $1
