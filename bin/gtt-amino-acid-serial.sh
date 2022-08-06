#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
ORANGE='\033[0;33m'
NC='\033[0m'

tmp_dir=$2
hmm_file=$3
amino_acid_genomes_total=$4
num_cpus=$5
hmm_target_genes_total=$6
output_dir=$7
best_hit_mode=$8
additional_pfam_targets=${9}
ko_targets=${10}
target_KOs=${11}

# looping through the lines of the provided [-f] file (this loop operates on one genome at a time)
while IFS=$'\t' read -r -a file
do

    ### kill backstop
    # if there is a problem on any iteration, exiting this subprocess and then exiting main script with report of problem assembly
    if [ -s ${tmp_dir}/kill_amino_acid_serial.problem ]; then
        exit
    fi


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
    echo $assembly >> ${tmp_dir}/amino_acid_genomes_list.tmp

    num=$((num+1)) # to track progress

    printf "   --------------------------------------------------------------------------   \n"
    printf "\tOn assembly ${GREEN}$assembly${NC}; Number $num of $amino_acid_genomes_total total.\n"
    printf "   --------------------------------------------------------------------------   \n\n"


    # filtering sequences by length to be sure none with > 99,999 amino acids are there, as this breaks hmmer (https://github.com/AstrobioMike/GToTree/issues/50 ; https://github.com/EddyRivasLab/hmmer/issues/244)
    gtt-filter-seqs-by-length -q -i ${file_location} -m 0 -M 99999 -o ${tmp_dir}/${assembly}_genes1.tmp

    ## renaming seqs to have assembly name (also to ensure simple headers)
    gtt-rename-fasta-headers -i ${tmp_dir}/${assembly}_genes1.tmp -w ${assembly} -o ${tmp_dir}/${assembly}_genes.tmp

    ## removing gunzipped genome file if it was gunzipped
    if [ $was_gzipped == "TRUE" ]; then
        rm -rf $file_location
    fi

    ## exiting here and reporting current input file if something is wrong with it and didn't get coding sequences
    if [ ! -s ${tmp_dir}/${assembly}_genes.tmp ]; then
        printf "$assembly" >> ${tmp_dir}/kill_amino_acid_serial.problem
        exit
    fi

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

    # adding NA for taxid so final table can still have the column and lineage for those that do have them
    taxid="NA"

    ## writing summary info to table ##
    printf "$assembly\t$file\t$taxid\t$num_SCG_hits\t$uniq_SCG_hits\t$perc_comp_rnd\t$perc_redund_rnd\n" >> ${output_dir}/Amino_acid_genomes_summary_info.tsv

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

    rm -rf ${tmp_dir}/${assembly}_genes1.tmp ${tmp_dir}/${assembly}_genes.tmp ${tmp_dir}/${assembly}_curr_hmm_hits.tmp
    rm -rf ${tmp_dir}/${assembly}_uniq_counts.tmp ${tmp_dir}/${assembly}_conservative_filtering_counts_tab.tmp
    rm -rf ${tmp_dir}/${assembly}_conservative_target_unique_hmm_names.tmp ${tmp_dir}/${assembly}_genes.tmp.ssi

done < $1
