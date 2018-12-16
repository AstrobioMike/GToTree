#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

tmp_dir=$(cat temp_dir_name.tmp)
len_cutoff=$(cat len_cutoff.tmp)

### filtering out sequences that are too long or too short ###

gtt-count-bases-per-seq -i ${tmp_dir}/${1}_hits.faa -o ${1}_Num_bps.tmp
cut -f2 ${1}_Num_bps.tmp > ${1}_lengths.tmp
median=$(gtt-get-median.sh ${1}_lengths.tmp)
buff=$(echo "$median * $len_cutoff" | bc)
min_len=$(echo "$median - $buff" | bc)
min_len_rnd=$(printf "%.0f\n" $min_len)
max_len=$(echo "$median + $buff" | bc)
max_len_rnd=$(printf "%.0f\n" $max_len)

gtt-filter-seqs-by-length -i ${tmp_dir}/${1}_hits.faa -m $min_len_rnd -M $max_len_rnd -o ${tmp_dir}/${1}_hits_filtered.tmp > ${1}_filter.out.tmp

cat <(printf "\n    Filtering ${GREEN}${1}${NC} sequences by length...\n") ${1}_filter.out.tmp

rm ${1}_Num_bps.tmp ${1}_lengths.tmp ${1}_filter.out.tmp
