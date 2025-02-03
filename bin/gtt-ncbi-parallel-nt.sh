#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
ORANGE='\033[0;33m'
NC='\033[0m'

tmp_dir="$2"
hmm_file="$3"
num_cpus="$4"
hmm_target_genes_total="$5"
output_dir="$6"
best_hit_mode="$7"
additional_pfam_targets="$8"
http_flag="$9"
ko_targets="${10}"
target_KOs="${11}"
debug_flag="${12}"

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

# attempting to download genome fasta
curl --silent --retry 10 -o ${tmp_dir}/${assembly}_genome.tmp.gz "${base_link}/${end_path}_genomic.fna.gz"

# if http, then it pulls down a file still, it just isn't gzipped
# if ftp, no file is pulled down
# so to cover both cases, just making this need to be present and gzipped
if $(file ${tmp_dir}/${assembly}_genome.tmp.gz | grep -q gzip); then

    gunzip -f ${tmp_dir}/${assembly}_genome.tmp.gz

    prodigal -c -q -i ${tmp_dir}/${assembly}_genome.tmp -a ${tmp_dir}/${assembly}_genes1.faa.tmp -d ${tmp_dir}/${assembly}_genes1.fa.tmp > /dev/null

    tr -d '*' < ${tmp_dir}/${assembly}_genes1.faa.tmp > ${tmp_dir}/${assembly}_genes2.faa.tmp

    ## renaming seqs to have assembly name to avoid problems with odd characters and how hmmer parses and such
    gtt-rename-fasta-headers -i ${tmp_dir}/${assembly}_genes2.faa.tmp -w $assembly -o ${tmp_dir}/${assembly}_genes3.faa.tmp
    gtt-rename-fasta-headers -i ${tmp_dir}/${assembly}_genes1.fa.tmp -w $assembly -o ${tmp_dir}/${assembly}_genes.fa.tmp

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

    ### counting how many genes in this genome
    gene_count=$(grep -c ">" ${tmp_dir}/${assembly}_genes3.faa.tmp)

    # filtering sequences by length to be sure none with > 99,999 amino acids are there, as this breaks hmmer (https://github.com/AstrobioMike/GToTree/issues/50 ; https://github.com/EddyRivasLab/hmmer/issues/244)
    gtt-filter-seqs-by-length -q -i ${tmp_dir}/${assembly}_genes3.faa.tmp -m 0 -M 99999 -o ${tmp_dir}/${assembly}_genes.faa.tmp
        # don't need to filter the nucleotide fasta, as we will only be looking for some found in the amino-acid fasta

    ### running hmm search ###
    hmmsearch --cut_ga --cpu $num_cpus --tblout ${tmp_dir}/${assembly}_curr_hmm_hits.tmp $hmm_file ${tmp_dir}/${assembly}_genes.faa.tmp > /dev/null

    ### calculating % completion and redundancy ###
    for SCG in $(cat ${tmp_dir}/uniq_hmm_names.tmp)
    do
        grep -w -c "$SCG" ${tmp_dir}/${assembly}_curr_hmm_hits.tmp
    done > ${tmp_dir}/${assembly}_uniq_counts.tmp

    ## making list here of only those present in exactly 1 copy, to get count of "unique" SCG-hits
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
    printf "$assembly\t$downloaded_accession\t$ass_name\t$taxid\t$org_name\t$infraspecific_name\t$version_status\t$asm_level\t$num_SCG_hits\t$uniq_SCG_hits\t$perc_comp_rnd\t$perc_redund_rnd\n" >> ${output_dir}/NCBI_genomes_summary_info.tsv

    ### Pulling out hits for this genome (nucleotide as specified by user) ###
    target_genes_suffix="_genes.fa.tmp"

    # indexing
    esl-sfetch --index ${tmp_dir}/${assembly}${target_genes_suffix} > /dev/null

    # looping through and pulling out each first hit (hmm results tab is sorted by e-value):
    # if best-hit mode is off, then only pulling genes that were identified in exactly 1 copy
    if [ $best_hit_mode  == "false" ]; then

        for SCG in $(cat ${tmp_dir}/${assembly}_conservative_target_unique_hmm_names.tmp)
        do
            grep -w "$SCG" ${tmp_dir}/${assembly}_curr_hmm_hits.tmp | awk '!x[$3]++' | cut -f1 -d " " | esl-sfetch -f ${tmp_dir}/${assembly}${target_genes_suffix} - | sed "s/>.*$/>$assembly/" | sed 's/^Usage.*$//' | sed 's/^To see.*$//' | sed '/^$/d' >> ${tmp_dir}/${SCG}_hits.fa
        done

    # if best-hit mode is on, taking best hit
    else

        for SCG in $(cat ${tmp_dir}/uniq_hmm_names.tmp)
        do
            grep -w "$SCG" ${tmp_dir}/${assembly}_curr_hmm_hits.tmp | awk '!x[$3]++' | cut -f1 -d " " | esl-sfetch -f ${tmp_dir}/${assembly}${target_genes_suffix} - | sed "s/>.*$/>$assembly/" | sed 's/^Usage.*$//' | sed 's/^To see.*$//' | sed '/^$/d' >> ${tmp_dir}/${SCG}_hits.fa
        done

    fi


    ## searching for additional targets if provided
    # getting count of genes if there are additional targets
    if [ $ko_targets == "true" ] || [ $additional_pfam_targets == "true" ]; then

        gene_count=$(grep -c ">" ${tmp_dir}/${assembly}_genes.faa.tmp)

    fi

    ## KOs
    if [ $ko_targets == "true" ]; then

        gtt-run-kofamscan.sh ${assembly} ${tmp_dir}/${assembly}_genes.faa.tmp ${gene_count} ${output_dir}/KO_search_results/target-KOs.tsv ${output_dir}/KO_search_results/target_KO_profiles/ ${tmp_dir} ${output_dir} ${target_KOs}

    fi

    ## Pfams
    if [ $additional_pfam_targets == "true" ]; then

        gtt-run-additional-pfam-search.sh ${assembly} ${tmp_dir}/${assembly}_genes.faa.tmp ${gene_count} ${num_cpus} ${tmp_dir} ${output_dir}

    fi

    if [ $debug_flag == "true" ]; then
        if [ -s ${tmp_dir}/${assembly}_genes2.faa.tmp ]; then
            mv ${tmp_dir}/${assembly}_genes2.faa.tmp ${tmp_dir}/ncbi-downloads/${assembly}_protein.faa
        fi
        if [ -s ${tmp_dir}/${assembly}_genes1.fa.tmp ]; then
            mv ${tmp_dir}/${assembly}_genes1.fa.tmp ${tmp_dir}/ncbi-downloads/${assembly}_cds.fa
        fi
        if [ -s ${tmp_dir}/${assembly}_genome.tmp ]; then
            mv ${tmp_dir}/${assembly}_genome.tmp ${tmp_dir}/ncbi-downloads/${assembly}_genomic.fna
        fi
    fi

    rm -rf ${tmp_dir}/${assembly}_genes*.tmp*
    rm -rf ${tmp_dir}/${assembly}_curr_hmm_hits.tmp ${tmp_dir}/${assembly}_uniq_counts.tmp
    rm -rf ${tmp_dir}/${assembly}_conservative_filtering_counts_tab.tmp
    rm -rf ${tmp_dir}/${assembly}_conservative_target_unique_hmm_names.tmp

else
    printf "     ${ORANGE}******************************* ${NC}NOTICE ${ORANGE}*******************************${NC}  \n"
    printf "\t  $assembly's genome not successfully downloaded :(\n\n"
    printf "\t    Reported in \"${output_dir}/run_files/NCBI_accessions_not_downloaded.txt\"\n"
    printf "     ${ORANGE}************************************************************************ ${NC}\n"
    rm -rf ${tmp_dir}/${assembly}_genome.tmp.gz

    sleep 2

    echo $assembly >> ${output_dir}/NCBI_accessions_not_downloaded.txt

fi
