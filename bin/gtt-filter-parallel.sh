#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
ORANGE='\033[0;33m'
NC='\033[0m'

tmp_dir=$2
len_cutoff=$3

### filtering out sequences that are too long or too short ###

gtt-count-bases-per-seq -i ${tmp_dir}/${1}_hits.faa -o ${tmp_dir}/${1}_Num_bps.tmp
cut -f2 ${tmp_dir}/${1}_Num_bps.tmp > ${tmp_dir}/${1}_lengths.tmp
median=$(gtt-get-median.sh ${tmp_dir}/${1}_lengths.tmp)
buff=$(echo "$median * $len_cutoff" | bc)
min_len=$(echo "$median - $buff" | bc)
min_len_rnd=$(printf "%.0f\n" $min_len)
max_len=$(echo "$median + $buff" | bc)
max_len_rnd=$(printf "%.0f\n" $max_len)

gtt-filter-seqs-by-length -i ${tmp_dir}/${1}_hits.faa -m $min_len_rnd -M $max_len_rnd -o ${tmp_dir}/${1}_hits_filtered.tmp > ${tmp_dir}/${1}_filter.out.tmp

cat <(printf "\n    Filtering ${GREEN}${1}${NC} sequences by length...\n") ${tmp_dir}/${1}_filter.out.tmp

rm ${tmp_dir}/${1}_Num_bps.tmp ${tmp_dir}/${1}_lengths.tmp ${tmp_dir}/${1}_filter.out.tmp
