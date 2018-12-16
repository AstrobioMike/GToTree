#!/usr/bin/env bash

tmp_dir=$(cat temp_dir_name.tmp)
hmm_file=$(cat hmm_file_path.tmp)
num_cpus=$(cat num_cpus.tmp)
hmm_target_genes_total=$(cat hmm_target_genes_total.tmp)

assembly=$(echo "$1" | cut -f 1)
downloaded_accession=$(echo "$1" | cut -f 2)

# storing links to download stuff in variables (if a refseq was returned as identical to what was searched, using that)
if [ ${downloaded_accession:0:3} == "GCA" ]; then
  base_link=$(echo "$1" | cut -f 8) # FtpPath_GenBank
else
  base_link=$(echo "$1" | cut -f 9) # FtpPath_RefSeq
fi

end_path=$(basename $base_link)

curl --silent --connect-timeout 10 --max-time 10 --retry 10 --retry-max-time 30 -o ${assembly}_report1.tmp "${base_link}/${end_path}_assembly_report.txt"
curl --silent --connect-timeout 10 --max-time 10 --retry 10 --retry-max-time 30 -o ${assembly}_genes.tmp.gz "${base_link}/${end_path}_protein.faa.gz"

if [ -s ${assembly}_report1.tmp ]; then
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
        else
          printf "     ${RED}******************************* ${NC}NOTICE ${RED}*******************************${NC}  \n"
          printf "\t  $assembly's genes nor genome downloaded properly :(\n\n"
          printf "\t    Reported in \"NCBI_accessions_not_downloaded.txt\"\n"
          printf "     ${RED}************************************************************************ ${NC}\n"
          rm -rf ${assembly}_report1.tmp ${assembly}_genes.tmp.gz
          sleep 3
          echo $assembly >> NCBI_accessions_not_downloaded.txt
        fi
    fi


# fixing the stupid carriage returns that for some reason NCBI assembly reports have in them...
    tr "\r" "\n" < ${assembly}_report1.tmp > ${assembly}_report.tmp

    # storing more info about the assembly to write out into ncbi-derived-genome summary file (for each setting to NA if not found)
    ass_name=$(grep "Assembly name:" ${assembly}_report.tmp | cut -f2 -d ":" | sed 's/^ *//')
    if [ -z "$ass_name" ]; then ass_name="NA"; fi
    org_name=$(grep "Organism name:" ${assembly}_report.tmp | cut -f2 -d ":" | sed 's/^ *//')
    if [ -z "$org_name" ]; then org_name="NA"; fi
    infraspecific_name=$(grep "Infraspecific name:" ${assembly}_report.tmp | cut -f2 -d ":" | sed 's/^ *//')
    if [ -z "$infraspecific_name" ]; then infraspecific_name="NA"; fi
    taxid=$(grep "Taxid:" ${assembly}_report.tmp | cut -f2 -d ":" | sed 's/^ *//')
    if [ -z "$taxid" ]; then taxid="NA"; fi

    printf "   ${GREEN}$assembly${NC}\n"
    printf "      Performing HMM search...\n"
      
    ### running hmm search ###
    # cp /Users/Mike_Lee/Documents/github/GToTree/hmm_sets/Hug_et_al.hmm ${assembly}.hmm
    hmmsearch --cut_ga --cpu $num_cpus --tblout ${assembly}_curr_hmm_hits.tmp $hmm_file ${assembly}_genes.tmp > /dev/null
    wait 
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
    printf "$assembly\t$downloaded_accession\t$ass_name\t$taxid\t$org_name\t$infraspecific_name\t$num_SCG_hits\t$perc_comp_rnd\t$perc_redund_rnd\n" >> NCBI_genomes_summary_info.tsv

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
    printf "\t  $assembly did not download properly :(\n\n"
    printf "\t    Reported in \"NCBI_accessions_not_downloaded.txt\"\n"
    printf "     ${RED}************************************************************************ ${NC}\n"
    rm -rf ${assembly}_report1.tmp ${assembly}_genes.tmp.gz
    sleep 3
    echo $assembly >> NCBI_accessions_not_downloaded.txt
fi

