#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
ORANGE='\033[0;33m'
NC='\033[0m'

tmp_dir=$2
hmm_file=$3
genbank_genomes_total=$4
num_cpus=$5
hmm_target_genes_total=$6
output_dir=$7
best_hit_mode=$8
additional_pfam_targets=$9
ko_targets=${10}
target_KOs=${11}

num=0

rm -rf ${output_dir}/run_files/Genbank_files_with_no_CDSs.txt # deleting if file exists

# looping through the lines of the provided [-g] file (this loop operates on one genome at a time)
while IFS=$'\t' read -r -a file

do

    ## checking if gzipped, gunzipping if so, and setting assembly name and file location variable either way
    if $(file $file | grep -q "gzip"); then
        was_gzipped=TRUE # setting variable to be able to check and remove gunzipped file afterwards
        file_location=${file%.*}
        gunzip -f -c $file > $file_location
        assembly="$(basename ${file_location%.*})"
    else
        file_location=$file
        assembly="$(basename ${file%.*})"
        was_gzipped=FALSE
    fi


    # adding assembly to ongoing genomes list
    echo $assembly >> ${tmp_dir}/genbank_genomes_list.tmp

    num=$((num+1)) # to track progress

    printf "   --------------------------------------------------------------------------   \n"
    printf "\tOn assembly ${GREEN}$assembly${NC}; Number $num of $genbank_genomes_total total.\n"
    printf "   --------------------------------------------------------------------------   \n\n"

    # storing more info about the assembly if it's present in the genbank file:
    if grep -q "ORGANISM" $file_location; then 
        org_name=$(grep -m1 "ORGANISM" $file_location | tr -s " " | cut -f3- -d " " | tr "[ ./\\]" "_" | tr -s "_")
    else
        org_name="NA"
    fi

    if grep -q "strain=" $file_location; then 
        strain=$(grep -m1 "strain=" $file_location | tr -s " " | cut -f 2 -d '"')
    else
        strain="NA"
    fi

    if grep -q "taxon" $file_location; then
        taxid=$(grep -m1 "taxon" $file_location | cut -f2 -d ":" | tr -d '"')
    else
        taxid="NA"
    fi

    # extracting AA coding sequences from genbank file
    gtt-genbank-to-AA-seqs -i $file_location -o ${tmp_dir}/${assembly}_genes2.tmp 2> /dev/null

    # checking that the file had CDS annotations, if not running prodigal
    if [ ! -s ${tmp_dir}/${assembly}_genes2.tmp ]; then

        printf "  ${ORANGE}********************************** ${NC}NOTICE ${ORANGE}**********************************${NC}  \n"
        printf "\t  This genbank file doesn't appear to have CDS annotations, so we are\n"
        printf "\t  identifying coding sequences with prodigal.\n\n"

        printf "\t    Reported in \"${output_dir}/run_files/Genbank_files_with_no_CDSs.txt\".\n"
        printf "  ${ORANGE}****************************************************************************${NC}  \n\n"

        echo "$file" >> ${output_dir}/run_files/Genbank_files_with_no_CDSs.txt
        rm -rf ${tmp_dir}/${assembly}_genes2.tmp

        # pulling out full nucleotide fasta from genbank file
        gtt-genbank-to-fasta -i $file_location -o ${tmp_dir}/${assembly}_fasta.tmp 2> /dev/null

        # running prodigal
        echo "prodigal used" > ${tmp_dir}/prodigal_used # marking so can add to citations list reported at end
        prodigal -c -q -i ${tmp_dir}/${assembly}_fasta.tmp -a ${tmp_dir}/${assembly}_genes1.tmp > /dev/null 2> ${file_location}_prodigal.stderr

        if grep -q "at least 100000 bases for training." ${file_location}_prodigal.stderr; then
            printf "$assembly\n" >> ${tmp_dir}/kill_genbank_serial.prodigal
            rm -rf ${file_location}_prodigal.stderr
            exit
        else
            rm -rf ${file_location}_prodigal.stderr
        fi

        tr -d '*' < ${tmp_dir}/${assembly}_genes1.tmp > ${tmp_dir}/${assembly}_genes2.tmp

    fi

    ## removing gunzipped genome file if it was gunzipped
    if [ $was_gzipped == "TRUE" ]; then
        rm -rf $file_location
    fi


    ## renaming seqs to have assembly name
    gtt-rename-fasta-headers -i ${tmp_dir}/${assembly}_genes2.tmp -w $assembly -o ${tmp_dir}/${assembly}_genes.tmp


    ### counting how many genes in this genome
    gene_count=$(grep -c ">" ${tmp_dir}/${assembly}_genes.tmp)

    printf "      Performing HMM search...\n"
      
    ### running hmm search ###
    hmmsearch --cut_ga --cpu $num_cpus --tblout ${tmp_dir}/${assembly}_curr_hmm_hits.tmp $hmm_file ${tmp_dir}/${assembly}_genes.tmp > /dev/null

    ### calculating % completion and redundancy ###
    for SCG in $(cat ${tmp_dir}/uniq_hmm_names.tmp)
    do
        grep -w -c "$SCG" ${tmp_dir}/${assembly}_curr_hmm_hits.tmp
    done > ${tmp_dir}/${assembly}_uniq_counts.tmp

    ## making list here of only those present in exactly 1 copy
    paste ${tmp_dir}/uniq_hmm_names.tmp ${tmp_dir}/${assembly}_uniq_counts.tmp > ${tmp_dir}/${assembly}_conservative_filtering_counts_tab.tmp
    awk -F "\t" ' $2 == 1 ' ${tmp_dir}/${assembly}_conservative_filtering_counts_tab.tmp | cut -f 1 > ${tmp_dir}/${assembly}_conservative_target_unique_hmm_names.tmp
    uniq_SCG_hits=$(wc -l ${tmp_dir}/${assembly}_conservative_target_unique_hmm_names.tmp | sed 's/^ *//' | cut -f 1 -d " ")


    ## adding SCG-hit counts to table
    paste <(printf $assembly) <(printf %s "$(cat ${tmp_dir}/${assembly}_uniq_counts.tmp | tr "\n" "\t" | sed 's/.$//')") >> ${output_dir}/SCG_hit_counts.tsv

    # total number of unique SCG hits
    num_SCG_hits=$(awk ' $1 > 0 ' ${tmp_dir}/${assembly}_uniq_counts.tmp | wc -l | tr -s " " | cut -f2 -d " ")
    num_SCG_redund=$(awk '{ if ($1 == 0) { print $1 } else { print $1 - 1 } }' ${tmp_dir}/${assembly}_uniq_counts.tmp | awk '{ sum += $1 } END { print sum }')

    perc_comp=$(echo "$num_SCG_hits / $hmm_target_genes_total * 100" | bc -l)
    perc_comp_rnd=$(printf "%.2f\n" $perc_comp)
    perc_redund=$(echo "$num_SCG_redund / $hmm_target_genes_total * 100" | bc -l)
    perc_redund_rnd=$(printf "%.2f\n" $perc_redund)

    # want to put a notice out if estimated redundancy is greater than 10
    # needs to be an integer for bash comparison, so multiplying by 100 first

    mult_perc_redund_rnd=$(echo "$perc_redund_rnd * 100" | bc | cut -f 1 -d ".")

    printf "        Found $num_SCG_hits of the targeted $hmm_target_genes_total genes.\n"

    if [ ${mult_perc_redund_rnd} -ge 1000 ]; then

        printf "        Est. %% comp: ${perc_comp_rnd}; Est. %% redund: ${RED}${perc_redund_rnd}${NC}\n\n"

        printf "  ${ORANGE}********************************** ${NC}NOTICE ${ORANGE}**********************************${NC}  \n"
        printf "   Estimated redundancy of this genome based on the specified HMMs is ${RED}${perc_redund_rnd}%%${NC}.\n"
        printf "   While there are no \"golden\" cutoff values for these things, typically\n"
        printf "   going over 10%% (if bacterial/archaeal) is getting into the questionable range.\n"
        printf "   You may want to consider taking a closer look and/or removing it from the\n"
        printf "   from the input genomes.\n\n"

        printf "   Reported in \"${output_dir}/run_files/Genomes_with_questionable_redund_estimates.tsv\".\n"
        printf "  ${ORANGE}****************************************************************************${NC}  \n\n"

        # writing to table of genomes with questionable redundancy estimates
        printf "$assembly\t$num_SCG_hits\t$uniq_SCG_hits\t$perc_comp_rnd\t$perc_redund_rnd\n" >> ${tmp_dir}/Genomes_with_questionable_redundancy_estimates.tmp

    else
        printf "        Est. %% comp: ${perc_comp_rnd}; Est. %% redund: ${perc_redund_rnd}\n\n"

    fi


    ## writing summary info to table ##
    printf "$assembly\t$file\t$taxid\t$org_name\t$strain\t$num_SCG_hits\t$uniq_SCG_hits\t$perc_comp_rnd\t$perc_redund_rnd\n" >> ${output_dir}/Genbank_genomes_summary_info.tsv

    ### Pulling out hits for this genome ###
    # looping through SCGs and pulling out each first hit (hmm results tab is sorted by e-value):
    esl-sfetch --index ${tmp_dir}/${assembly}_genes.tmp > /dev/null
        
    # if best-hit mode is off, then only pulling genes that were identified in exactly 1 copy
    if [ $best_hit_mode  == "false" ]; then

        for SCG in $(cat ${tmp_dir}/${assembly}_conservative_target_unique_hmm_names.tmp)
        do
            grep -w "$SCG" ${tmp_dir}/${assembly}_curr_hmm_hits.tmp | awk '!x[$3]++' | cut -f1 -d " " | esl-sfetch -f ${tmp_dir}/${assembly}_genes.tmp - | sed "s/>.*$/>$assembly/" | sed 's/^Usage.*$//' | sed 's/^To see.*$//' | sed '/^$/d' >> ${tmp_dir}/${SCG}_hits.faa
        done

    # if best-hit mode is on, taking best hit
    else

        for SCG in $(cat ${tmp_dir}/uniq_hmm_names.tmp)
        do
            grep -w "$SCG" ${tmp_dir}/${assembly}_curr_hmm_hits.tmp | awk '!x[$3]++' | cut -f1 -d " " | esl-sfetch -f ${tmp_dir}/${assembly}_genes.tmp - | sed "s/>.*$/>$assembly/" | sed 's/^Usage.*$//' | sed 's/^To see.*$//' | sed '/^$/d' >> ${tmp_dir}/${SCG}_hits.faa
        done

    fi


    ## searching for additional targets if provided
    # getting count of genes if there are additional targets
    if [ $ko_targets == "true" ] || [ $additional_pfam_targets == "true" ]; then

        gene_count=$(grep -c ">" ${tmp_dir}/${assembly}_genes.tmp)

    fi

    ## KOs
    if [ $ko_targets == "true" ]; then

        gtt-run-kofamscan.sh ${assembly} ${tmp_dir}/${assembly}_genes.tmp ${gene_count} ${output_dir}/KO_search_results/target-KOs.tsv ${output_dir}/KO_search_results/target_KO_profiles/ ${tmp_dir} ${output_dir} ${target_KOs}

    fi

    ## Pfams
    if [ $additional_pfam_targets == "true" ]; then

        gtt-run-additional-pfam-search.sh ${assembly} ${tmp_dir}/${assembly}_genes.tmp ${gene_count} ${num_cpus} ${tmp_dir} ${output_dir}

    fi

    rm -rf ${tmp_dir}/${assembly}_*.tmp ${tmp_dir}/${assembly}_genes.tmp.ssi

done < $1
