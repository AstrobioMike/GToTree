#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
ORANGE='\033[0;33m'
NC='\033[0m'

tmp_dir=$2
hmm_file=$3
num_cpus=$4
hmm_target_genes_total=$5
output_dir=$6
best_hit_mode=$7
additional_pfam_targets=$8
http_flag=$9

assembly=$(echo "$1" | cut -f 1)
downloaded_accession=$(echo "$1" | cut -f 2)

# storing and building links
if [ "$http_flag" == 'false' ]; then
    base_link=$(echo "$1" | cut -f 9)
else
    base_link=$(echo "$1" | cut -f 9 | sed 's/^ftp/https/')
fi

# checking link was actually present (sometimes, very rarely, it is not there)
# if not there, attempting to build ourselves
if [ $base_link == "na" ] || [ -z $base_link ]; then

    if [ "$http_flag" == 'false' ]; then
        p1=$(printf "ftp://ftp.ncbi.nlm.nih.gov/genomes/all")
    else
        p1=$(printf "https://ftp.ncbi.nlm.nih.gov/genomes/all")
    fi

    # checking if GCF or GCA
    if [[ $assembly == "GCF"* ]]; then 
        p2="GCF"
    else
        p2="GCA"
    fi
    
    p3=$(echo $assembly | cut -f 2 -d "_" | cut -c 1-3)
    p4=$(echo $assembly | cut -f 2 -d "_" | cut -c 4-6)
    p5=$(echo $assembly | cut -f 2 -d "_" | cut -c 7-9)

    ass_name=$(echo "$1" | cut -f 3)
    end_path=$(paste -d "_" <(echo "$assembly") <(echo "$ass_name"))

    base_link=$(paste -d "/" <(echo "$p1") <(echo "$p2") <(echo "$p3") <(echo "$p4") <(echo "$p5") <(echo "$end_path"))

else

    end_path=$(basename $base_link)

fi

printf "   --------------------------------------------------------------------------   \n\n"
printf "     Genome: ${GREEN}$assembly${NC}\n"

curl --silent --retry 10 -o ${tmp_dir}/${assembly}_genes2.tmp.gz "${base_link}/${end_path}_protein.faa.gz"

if [ -s ${tmp_dir}/${assembly}_genes2.tmp.gz ]; then
    gunzip ${tmp_dir}/${assembly}_genes2.tmp.gz
    # renaming headers to avoid problems with odd characters and how hmmer parses and such
    gtt-rename-fasta-headers -i ${tmp_dir}/${assembly}_genes2.tmp -w $assembly -o ${tmp_dir}/${assembly}_genes.tmp

else # trying to get assembly if there were no gene annotations available
    curl --silent --retry 10 -o ${tmp_dir}/${assembly}_genome.tmp.gz "${base_link}/${end_path}_genomic.fna.gz"

    if [ -s ${tmp_dir}/${assembly}_genome.tmp.gz ]; then

      gunzip ${tmp_dir}/${assembly}_genome.tmp.gz

      printf "  ${ORANGE}********************************** ${NC}NOTICE ${ORANGE}**********************************${NC}  \n"
      printf "   $assembly doesn't appear to have gene annotations.\n\n"
      printf "   Downloaded the genome and identifying CDSs with prodigal.\n"
      printf "  ${ORANGE}****************************************************************************${NC}  \n\n"

      printf "      Getting coding seqs...\n\n"
      echo "prodigal used" > ${tmp_dir}/prodigal_used # marking so can add to citations list reported at end
      prodigal -c -q -i ${tmp_dir}/${assembly}_genome.tmp -a ${tmp_dir}/${assembly}_genes1.tmp > /dev/null
      tr -d '*' < ${tmp_dir}/${assembly}_genes1.tmp > ${tmp_dir}/${assembly}_genes2.tmp

      ## renaming seqs to have assembly name
      gtt-rename-fasta-headers -i ${tmp_dir}/${assembly}_genes2.tmp -w $assembly -o ${tmp_dir}/${assembly}_genes.tmp
    fi
fi

if [ -s ${tmp_dir}/${assembly}_genes.tmp ]; then

    # storing more info about the assembly to write out into ncbi-derived-genome summary file (for each setting to NA if not found)
    ass_name=$(echo "$1" | cut -f 3)
    if [ -z "$ass_name" ]; then ass_name="NA"; fi
    org_name=$(echo "$1" | cut -f 5)
    if [ -z "$org_name" ]; then org_name="NA"; fi
    infraspecific_name=$(echo "$1" | cut -f 6)
    if [ -z "$infraspecific_name" ]; then infraspecific_name="NA"; fi
    taxid=$(echo "$1" | cut -f 4)
    if [ -z "$taxid" ]; then taxid="NA"; fi
    version_status=$(echo "$1" | cut -f 7)
    if [ -z "$version_status" ]; then version_status="NA"; fi
    asm_level=$(echo "$1" | cut -f 8)
    if [ -z "$asm_level" ]; then asm_level="NA"; fi



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
    paste <(printf $assembly) <(printf %s "$(cat ${tmp_dir}/${assembly}_uniq_counts.tmp | tr "\n" "\t")") >> ${output_dir}/SCG_hit_counts.tsv

    num_SCG_hits=$(awk ' $1 > 0 ' ${tmp_dir}/${assembly}_uniq_counts.tmp | wc -l | tr -s " " | cut -f2 -d " ")
    num_SCG_redund=$(awk '{ if ($1 == 0) { print $1 } else { print $1 - 1 } }' ${tmp_dir}/${assembly}_uniq_counts.tmp | awk '{ sum += $1 } END { print sum }')

    perc_comp=$(echo "$uniq_SCG_hits / $hmm_target_genes_total * 100" | bc -l)
    perc_comp_rnd=$(printf "%.2f\n" $perc_comp)
    perc_redund=$(echo "$num_SCG_redund / $hmm_target_genes_total * 100" | bc -l)
    perc_redund_rnd=$(printf "%.2f\n" $perc_redund)

    ### want to put an explicit notice out if estimated redundancy is greater than 10%
    # needs to be an integer for bash comparison, so multiplying by 100 first

    mult_perc_redund_rnd=$(echo "$perc_redund_rnd * 100" | bc | cut -f 1 -d ".")

    printf "             Found $num_SCG_hits of the targeted $hmm_target_genes_total genes.\n"

    if [ ${mult_perc_redund_rnd} -ge 1000 ]; then
        printf "             Est. %% comp: ${perc_comp_rnd}; Est. %% redund: ${RED}${perc_redund_rnd}${NC}\n\n"


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
        printf "             Est. %% comp: ${perc_comp_rnd}; Est. %% redund: ${perc_redund_rnd}\n\n"
    fi

    ## writing summary info to table ##
    printf "$assembly\t$downloaded_accession\t$ass_name\t$taxid\t$org_name\t$infraspecific_name\t$version_status\t$asm_level\t$num_SCG_hits\t$uniq_SCG_hits\t$perc_comp_rnd\t$perc_redund_rnd\n" >> ${output_dir}/NCBI_genomes_summary_info.tsv

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
        paste <( printf "${assembly}\t${downloaded_accession}\t${gene_count}" ) <(printf %s "$(cat ${tmp_dir}/${assembly}_hit_counts.tmp | tr "\n" "\t") " )  >> ${output_dir}/additional_pfam_search_results/Additional_Pfam_hit_counts.tsv

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


else
    printf "  ${ORANGE}********************************** ${NC}NOTICE ${ORANGE}**********************************${NC}  \n"
    printf "   $assembly did not download properly :(\n\n"
    printf "   Reported in \"${output_dir}/run_files/NCBI_accessions_not_downloaded.txt\"\n"
    printf "  ${ORANGE}****************************************************************************${NC}  \n\n"
    rm -rf ${tmp_dir}/${assembly}_report1.tmp ${tmp_dir}/${assembly}_genes.tmp.gz
    sleep 3
    echo $assembly >> ${output_dir}/NCBI_accessions_not_downloaded.txt
fi
