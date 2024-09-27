#!/usr/bin/env bash

tmp_dir=${1}
output_dir=${2}

# base_link="https://pfam.xfam.org/family/"
    # base link updated Oct-2022 when pfam hosting shifted to interpro
base_link="https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/"

# starting table of what was requested and what was found (version-wise, the latest is always pulled from Pfam when downloading as below)
printf "requested_Pfam\tpulled_Pfam\n" > ${output_dir}/Pfam_search_results/info/Requested-and-pulled.tsv

for initial_target in $(cat ${tmp_dir}/uniq_pfam_targets.tmp)
do


    # getting target without version specified if there was one (including a version doesn't work anymore since pfam-hosting shifted to interpro)
    target=$(echo ${initial_target} | cut -f 1 -d ".")

    # --insecure flag added on 29-Nov-2020, due to pfam certificate being invalid (https://github.com/AstrobioMike/GToTree/issues/28)
    curl --insecure --silent --retry 10 -o ${tmp_dir}/${target}.hmm.gz "${base_link}${target}?annotation=hmm"
    gunzip ${tmp_dir}/${target}.hmm.gz

    if [ -s ${tmp_dir}/${target}.hmm ]; then
        # getting accession pulled (to account for current version on Pfam as compared to what was searched)
        actual_target=$(grep -m1 "^ACC" ${tmp_dir}/${target}.hmm | tr -s " " "\t" | cut -f 2)
        printf "$actual_target\n" >> ${tmp_dir}/actual_pfam_targets.tmp

        if [ $initial_target != $actual_target ]; then
            mv ${tmp_dir}/${target}.hmm ${tmp_dir}/${actual_target}.hmm
        fi

        cat ${tmp_dir}/${actual_target}.hmm >> ${tmp_dir}/all_pfam_targets.hmm

        # adding searched and pulled to info table (meaning which versions of a Pfam)
        printf "${initial_target}\t${actual_target}\n" >> ${output_dir}/Pfam_search_results/info/Requested-and-pulled.tsv

    else # aborting if any of the pfam targets couldn't be pulled successfully
        printf "\n  ${RED}One of the target Pfams could not be successfully downloaded :(${NC}\n"
        printf "\n  The problem child was ${target}.\n\n"
        printf "\nExiting for now.\n\n"

        rm -rf ${output_dir}
        # removing temp directory unless debug mode on
        if [ $debug_flag == 'false' ]; then
             rm -rf $tmp_dir
        fi

        exit

    fi

done

# starting the main results table which will have the following as its header:
paste <(printf "assembly_id\ttotal_gene_count") <(printf %s "$(cat ${tmp_dir}/actual_pfam_targets.tmp | tr "\n" "\t")") > ${output_dir}/Pfam_search_results/Pfam-hit-counts.tsv

# copying over the additional pfam-target hmms
cp ${tmp_dir}/all_pfam_targets.hmm ${output_dir}/Pfam_search_results/target_Pfam_profiles/all-additional-pfam-targets.hmm
