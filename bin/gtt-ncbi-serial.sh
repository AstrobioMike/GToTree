#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

tmp_dir=$2
hmm_file=$3
NCBI_remaining_genomes_total=$4
num_cpus=$5
hmm_target_genes_total=$6
output_dir=$7
best_hit_mode=$8

num=0

while IFS=$'\t' read -r -a curr_line

do

    assembly="${curr_line[0]}"
    downloaded_accession="${curr_line[1]}"
    num=$((num+1))

    printf "   --------------------------------------------------------------------------   \n"
    printf "\tOn assembly ${GREEN}$assembly${NC}; Number $num of $NCBI_remaining_genomes_total total.\n"
    printf "   --------------------------------------------------------------------------   \n\n"

    # storing and building links
    base_link="${curr_line[8]}"
    end_path=$(basename $base_link)

    # attempting to download genes for assembly
    curl --silent --retry 10 -o ${tmp_dir}/${assembly}_genes2.tmp.gz "${base_link}/${end_path}_protein.faa.gz"

    if [ -s ${tmp_dir}/${assembly}_genes2.tmp.gz ]; then
        gunzip ${tmp_dir}/${assembly}_genes2.tmp.gz
        # renaming headers to avoid problems with odd characters and how hmmer parses and such
        gtt-rename-fasta-headers -i ${tmp_dir}/${assembly}_genes2.tmp -w $assembly -o ${tmp_dir}/${assembly}_genes.tmp

    else # trying to get assembly if there were no gene annotations available
        curl --silent --retry 10 -o ${tmp_dir}/${assembly}_genome.tmp.gz "${base_link}/${end_path}_genomic.fna.gz"

        if [ -s ${tmp_dir}/${assembly}_genome.tmp.gz ]; then

            gunzip ${tmp_dir}/${assembly}_genome.tmp.gz

            printf "  ${RED}********************************** ${NC}NOTICE ${RED}**********************************${NC}  \n"
            printf "   $assembly doesn't appear to have gene annotations.\n\n"
            printf "   Downloaded the genome and identifying CDSs with prodigal.\n"
            printf "  ${RED}****************************************************************************${NC}  \n\n"

            printf "      Getting coding seqs...\n\n"
            prodigal -c -q -i ${tmp_dir}/${assembly}_genome.tmp -a ${tmp_dir}/${assembly}_genes1.tmp > /dev/null
            tr -d '*' < ${tmp_dir}/${assembly}_genes1.tmp > ${tmp_dir}/${assembly}_genes2.tmp

            ## renaming seqs to have assembly name
            gtt-rename-fasta-headers -i ${tmp_dir}/${assembly}_genes2.tmp -w $assembly -o ${tmp_dir}/${assembly}_genes.tmp
        fi
    fi

    if [ -s ${tmp_dir}/${assembly}_genes.tmp ]; then

        # storing more info about the assembly to write out into ncbi-derived-genome summary file (for each setting to NA if not found)
        ass_name="${curr_line[2]}"
        if [ -z "$ass_name" ]; then ass_name="NA"; fi
        org_name="${curr_line[4]}"
        if [ -z "$org_name" ]; then org_name="NA"; fi
        infraspecific_name="${curr_line[5]}"
        if [ -z "$infraspecific_name" ]; then infraspecific_name="NA"; fi
        taxid="${curr_line[3]}"
        if [ -z "$taxid" ]; then taxid="NA"; fi
        version_status="${curr_line[6]}"
        if [ -z "$version_status" ]; then version_status="NA"; fi
        asm_level="${curr_line[7]}"
        if [ -z "$asm_level" ]; then asm_level="NA"; fi

        printf "      Performing HMM search...\n"

        ### running hmm search ###
        hmmsearch --cut_ga --cpu $num_cpus --tblout ${tmp_dir}/${assembly}_curr_hmm_hits.tmp $hmm_file ${tmp_dir}/${assembly}_genes.tmp > /dev/null

        ### calculating % completion and redundancy ###
        for SCG in $(cat ${tmp_dir}/uniq_hmm_names.tmp)
        do
            grep -w -c "$SCG" ${tmp_dir}/${assembly}_curr_hmm_hits.tmp
        done > ${tmp_dir}/${assembly}_uniq_counts.tmp

        ## making list here of only those present in exactly 1 copy, unless running in "best-hit mode" ("-B" flag specified)
        if [ $best_hit_mode == "false" ]; then
            paste ${tmp_dir}/uniq_hmm_names.tmp ${tmp_dir}/${assembly}_uniq_counts.tmp > ${tmp_dir}/${assembly}_conservative_filtering_counts_tab.tmp
            awk -F "\t" ' $2 == 1 ' ${tmp_dir}/${assembly}_conservative_filtering_counts_tab.tmp | cut -f 1 > ${tmp_dir}/${assembly}_conservative_target_unique_hmm_names.tmp
        fi

        ## adding SCG-hit counts to table
        paste <(printf $assembly) <(printf %s "$(cat ${tmp_dir}/${assembly}_uniq_counts.tmp | tr "\n" "\t")") >> ${output_dir}/All_genomes_SCG_hit_counts.tsv

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

            printf "  ${RED}********************************** ${NC}NOTICE ${RED}**********************************${NC}  \n"
            printf "   Estimated redundancy of this genome based on the specified HMMs is ${RED}${perc_redund_rnd}%%${NC}.\n"
            printf "   While there are no \"golden\" cutoff values for these things, typically\n"
            printf "   going over 10%% is getting into the questionable range. You may want to\n"
            printf "   consider taking a closer look and/or removing it from the input genomes.\n\n"

            printf "   Reported in \"${output_dir}/Genomes_with_questionable_redund_estimates.tsv\".\n"
            printf "  ${RED}****************************************************************************${NC}  \n\n"

            # writing to table of genomes with questionable redundancy estimates
            printf "$assembly\t$num_SCG_hits\t$perc_comp_rnd\t$perc_redund_rnd\n" >> ${tmp_dir}/Genomes_with_questionable_redundancy_estimates.tmp

        else
            printf "        Est. %% comp: ${perc_comp_rnd}; Est. %% redund: ${perc_redund_rnd}\n\n"

        fi


        ## writing summary info to table ##
        printf "$assembly\t$downloaded_accession\t$ass_name\t$taxid\t$org_name\t$infraspecific_name\t$version_status\t$asm_level\t$num_SCG_hits\t$perc_comp_rnd\t$perc_redund_rnd\n" >> ${output_dir}/NCBI_genomes_summary_info.tsv

        ### Pulling out hits for this genome ###
        # looping through and pulling out each first hit (hmm results tab is sorted by e-value):

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

        rm -rf ${tmp_dir}/${assembly}_*.tmp ${tmp_dir}/${assembly}_genes.tmp.ssi

    else
        printf "     ${RED}******************************* ${NC}NOTICE ${RED}*******************************${NC}  \n"
        printf "\t  $assembly's genes nor genome downloaded properly :(\n\n"
        printf "\t    Reported in \"${output_dir}/NCBI_accessions_not_downloaded.txt\"\n"
        printf "     ${RED}************************************************************************ ${NC}\n"
        rm -rf ${tmp_dir}/${assembly}_report1.tmp ${tmp_dir}/${assembly}_genes.tmp.gz
        sleep 3
        echo $assembly >> ${output_dir}/NCBI_accessions_not_downloaded.txt

    fi

done < $1
