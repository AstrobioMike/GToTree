#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

tmp_dir=$(cat temp_dir_name.tmp)
hmm_file=$(cat hmm_file_path.tmp)
NCBI_remaining_genomes_total=$(cat remaining_total_genomes.tmp)
num_cpus=$(cat num_cpus.tmp)
hmm_target_genes_total=$(cat hmm_target_genes_total.tmp)

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
    curl --silent --connect-timeout 10 --max-time 10 --retry 10 --retry-max-time 30 -o ${assembly}_genes.tmp.gz "${base_link}/${end_path}_protein.faa.gz"

    if [ -s ${assembly}_genes.tmp.gz ]; then
        gunzip ${assembly}_genes.tmp.gz
    else # trying to get assembly if there were no gene annotations available
        curl --silent --connect-timeout 10 --max-time 10 --retry 10 --retry-max-time 30 -o ${assembly}_genome.tmp.gz "${base_link}/${end_path}_genomic.fna.gz"
      
        if [ -s ${assembly}_genome.tmp.gz ]; then

            gunzip ${assembly}_genome.tmp.gz

            printf "     ${RED}******************************* ${NC}NOTICE ${RED}*******************************${NC}  \n"
            printf "\t  $assembly doesn't appear to have gene annotations.\n\n"
            printf "\t  Downloaded the genome and identifying CDSs with prodigal.\n"
            printf "     ${RED}********************************************************************** ${NC}\n\n"

            printf "      Getting coding seqs...\n\n"
            prodigal -c -q -i ${assembly}_genome.tmp -a ${assembly}_genes1.tmp > /dev/null
            tr -d '*' < ${assembly}_genes1.tmp > ${assembly}_genes2.tmp

            ## renaming seqs to have assembly name
            gtt-rename-fasta-headers -i ${assembly}_genes2.tmp -w $assembly -o ${assembly}_genes.tmp
        fi
    fi

    if [ -s ${assembly}_genes.tmp ]; then

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

        hmmsearch --cut_ga --cpu $num_cpus --tblout ${assembly}_curr_hmm_hits.tmp $hmm_file ${assembly}_genes.tmp > /dev/null

        ### calculating % completion and redundancy ###
        for SCG in $(cat ${tmp_dir}/uniq_hmm_names.tmp)
        do
            grep -w -c "$SCG" ${assembly}_curr_hmm_hits.tmp
        done > ${assembly}_uniq_counts.tmp

        num_SCG_hits=$(awk ' $1 > 0 ' ${assembly}_uniq_counts.tmp | wc -l | tr -s " " | cut -f2 -d " ")

        printf "        Found $num_SCG_hits of the targeted $hmm_target_genes_total SCGs.\n\n"

        num_SCG_redund=$(awk '{ if ($1 == 0) { print $1 } else { print $1 - 1 } }' ${assembly}_uniq_counts.tmp | awk '{ sum += $1 } END { print sum }')

        perc_comp=$(echo "$num_SCG_hits / $hmm_target_genes_total * 100" | bc -l)
        perc_comp_rnd=$(printf "%.2f\n" $perc_comp)
        perc_redund=$(echo "$num_SCG_redund / $hmm_target_genes_total * 100" | bc -l)
        perc_redund_rnd=$(printf "%.2f\n" $perc_redund)

        ## writing summary info to table ##
        printf "$assembly\t$downloaded_accession\t$ass_name\t$taxid\t$org_name\t$infraspecific_name\t$version_status\t$asm_level\t$num_SCG_hits\t$perc_comp_rnd\t$perc_redund_rnd\n" >> NCBI_genomes_summary_info.tsv

        ### Pulling out hits for this genome ###
        # making fasta file searchable to pull out the hits (Easel 0.45h June 2018)
        esl-sfetch --index ${assembly}_genes.tmp > /dev/null

        # looping through ribosomal proteins and pulling out each first hit (hmm results tab is sorted by e-value):
        # does as a separate loop just for clarity (well, in hopes of clarity)
        for SCG in $(cat ${tmp_dir}/uniq_hmm_names.tmp)
        do
            grep -w "$SCG" ${assembly}_curr_hmm_hits.tmp | awk '!x[$3]++' | cut -f1 -d " " | esl-sfetch -f ${assembly}_genes.tmp - | sed "s/>.*$/>$assembly/" | sed 's/^Usage.*$//' | sed 's/^To see.*$//' | sed '/^$/d' >> ${tmp_dir}/${SCG}_hits.faa
        done

        rm ${assembly}_*.tmp
        rm ${assembly}_*.tmp.ssi

    else
        printf "     ${RED}******************************* ${NC}NOTICE ${RED}*******************************${NC}  \n"
        printf "\t  $assembly's genes nor genome downloaded properly :(\n\n"
        printf "\t    Reported in \"NCBI_accessions_not_downloaded.txt\"\n"
        printf "     ${RED}************************************************************************ ${NC}\n"
        rm -rf ${assembly}_report1.tmp ${assembly}_genes.tmp.gz
        sleep 3
        echo $assembly >> NCBI_accessions_not_downloaded.txt

    fi

done < $1
