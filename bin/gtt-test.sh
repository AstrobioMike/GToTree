#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

printf "\n"

printf "  ${GREEN}Test data directory is located here:\n${NC}"
printf "    $TEST_DATA_DIR\n\n"

sleep 1

printf "  ${GREEN}Running test as:\n${NC}"
printf "    GToTree -a ${TEST_DATA_DIR}/ncbi_accessions.txt "'\\ \n'
printf "            -g ${TEST_DATA_DIR}/genbank_files.txt "'\\ \n'
printf "            -f ${TEST_DATA_DIR}/fasta_files.txt "'\\ \n'
printf "            -A ${TEST_DATA_DIR}/amino_acid_files.txt "'\\ \n'
printf "            -m ${TEST_DATA_DIR}/genome_to_id_map.tsv "'\\ \n'
printf "            -p ${TEST_DATA_DIR}/pfam_targets.txt "'\\ \n'
printf "            -H Universal -t -D -j 4 -o GToTree_test\n\n"

sleep 2

GToTree -a ${TEST_DATA_DIR}/ncbi_accessions.txt -g ${TEST_DATA_DIR}/genbank_files.txt -f ${TEST_DATA_DIR}/fasta_files.txt -A ${TEST_DATA_DIR}/amino_acid_files.txt -H Universal -m ${TEST_DATA_DIR}/genome_to_id_map.tsv -p ${TEST_DATA_DIR}/pfam_targets.txt -t -D -j 4 -o GToTree_test
