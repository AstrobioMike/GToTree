#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
ORANGE='\033[0;33m'
NC='\033[0m'

tmp_dir=$2
hmm_file=$3
fasta_genomes_total=$4
num_cpus=$5
hmm_target_genes_total=$6
output_dir=$7
best_hit_mode=$8
additional_pfam_targets=$9

# looping through the lines of the provided [-f] file (this loop operates on one genome at a time)
while IFS=$'\t' read -r -a file
do

    ### kill backstop
    # if there is a problem on any iteration, exiting this subprocess and then exiting main script with report of problem assembly
    if [ -s ${tmp_dir}/kill_fasta_serial.prodigal ]; then
        exit
    fi


    ## checking if gzipped, gunzipping if so, and setting assembly name and file location variable either way
    if $(file $file | grep -q "gzip"); then
        was_gzipped=TRUE # setting variable to be able to check and remove gunzipped file afterwards
        file_location=${file%.*}
        gunzip -c $file > $file_location
        assembly="$(basename ${file_location%.*})"
    else
        file_location=$file
        assembly="$(basename ${file%.*})"
        was_gzipped=FALSE
    fi

    # adding assembly to ongoing genomes list
    echo $assembly >> ${tmp_dir}/fasta_genomes_list.tmp

    num=$((num+1)) # to track progress

    printf "   --------------------------------------------------------------------------   \n"
    printf "\tOn assembly ${GREEN}$assembly${NC}; Number $num of $fasta_genomes_total total.\n"
    printf "   --------------------------------------------------------------------------   \n\n"

    printf "      Getting coding seqs...\n\n"

    ## running prodigal to get coding sequences
    prodigal -c -q -i $file_location -a ${tmp_dir}/${assembly}_genes1.tmp > /dev/null 2> ${file_location}_prodigal.stderr

    if grep -q "at least 100000 bases for training." ${file_location}_prodigal.stderr; then
        printf "$assembly\n" >> ${tmp_dir}/kill_fasta_serial.prodigal
        rm -rf ${file_location}_prodigal.stderr

        exit
    else
        rm -rf ${file_location}_prodigal.stderr
    fi

    tr -d '*' < ${tmp_dir}/${assembly}_genes1.tmp > ${tmp_dir}/${assembly}_genes2.tmp

    ## removing gunzipped genome file if it was gunzipped
    if [ $was_gzipped == "TRUE" ]; then
        rm -rf $file_location
    fi

    ## renaming seqs to have assembly name
    gtt-rename-fasta-headers -i ${tmp_dir}/${assembly}_genes2.tmp -w $assembly -o ${tmp_dir}/${assembly}_genes.tmp

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

    num_SCG_hits=$(awk ' $1 > 0 ' ${tmp_dir}/${assembly}_uniq_counts.tmp | wc -l | tr -s " " | cut -f2 -d " ")
    num_SCG_redund=$(awk '{ if ($1 == 0) { print $1 } else { print $1 - 1 } }' ${tmp_dir}/${assembly}_uniq_counts.tmp | awk '{ sum += $1 } END { print sum }')

    perc_comp=$(echo "$uniq_SCG_hits / $hmm_target_genes_total * 100" | bc -l)
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
        printf "   going over 10%% is getting into the questionable range. You may want to\n"
        printf "   consider taking a closer look and/or removing it from the input genomes.\n\n"

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
    printf "$assembly\t$file\t$taxid\t$num_SCG_hits\t$uniq_SCG_hits\t$perc_comp_rnd\t$perc_redund_rnd\n" >> ${output_dir}/Fasta_genomes_summary_info.tsv


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
    if [ $additional_pfam_targets == "true" ]; then

        ### counting how many genes in this genome
        gene_count=$(grep -c ">" ${tmp_dir}/${assembly}_genes.tmp)

        hmmsearch --cut_ga --cpu $num_cpus --tblout ${tmp_dir}/${assembly}_curr_hmm_hits.tmp ${tmp_dir}/all_targets.hmm ${tmp_dir}/${assembly}_genes.tmp > /dev/null

        ### getting counts of each target in this genome
        for target in $(cat ${tmp_dir}/actual_pfam_targets.tmp)
        do
            grep -w ${target} ${tmp_dir}/${assembly}_curr_hmm_hits.tmp | wc -l | sed 's/^ *//' >> ${tmp_dir}/${assembly}_hit_counts.tmp
        done

        ### writing results to main output file
        paste <( printf "${assembly}\tNA\t${gene_count}" ) <(printf %s "$(cat ${tmp_dir}/${assembly}_hit_counts.tmp | tr "\n" "\t") " )  >> ${output_dir}/additional_pfam_search_results/Additional_Pfam_hit_counts.tsv

        ### Pulling out hits to additional pfam targets for this genome ###
        for target in $(cat ${tmp_dir}/actual_pfam_targets.tmp)
        do
            if grep -w -q "$target" ${tmp_dir}/${assembly}_curr_hmm_hits.tmp; then

                grep -w "$target" ${tmp_dir}/${assembly}_curr_hmm_hits.tmp | cut -f 1 -d " " >> ${tmp_dir}/${assembly}_${target}_genes_of_int.tmp

                for gene in $(cat ${tmp_dir}/${assembly}_${target}_genes_of_int.tmp)
                do
                    echo $gene | esl-sfetch -f ${tmp_dir}/${assembly}_genes.tmp -
                done >> ${tmp_dir}/${assembly}_${target}_genes1.tmp

                gtt-append-fasta-headers -i ${tmp_dir}/${assembly}_${target}_genes1.tmp -w ${assembly}_${target} -o ${tmp_dir}/${assembly}_${target}_genes.tmp
            
                # adding to fasta of that target holding all genomes
                cat ${tmp_dir}/${assembly}_${target}_genes.tmp >> ${output_dir}/additional_pfam_search_results/${target}_hits.faa
            fi

        done
    fi

    rm -rf ${tmp_dir}/${assembly}_*.tmp ${tmp_dir}/${assembly}_genes.tmp.ssi

done < $1
