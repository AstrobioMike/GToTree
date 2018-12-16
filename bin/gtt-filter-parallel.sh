#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

tmp_dir=$(cat temp_dir_name.tmp)
len_cutoff=$(cat len_cutoff.tmp)

### filtering out sequences that are too long or too short ###
# creating short R scripts for getting filtering lengths
echo "lengths <- scan(\"${1}_lengths.tmp\", quiet=TRUE)
med <- median(lengths)
buff <- med*${len_cutoff}
minimum <- round(med - buff)
cat(minimum)" > ${tmp_dir}/${1}_min.R

echo "lengths <- scan(\"${1}_lengths.tmp\", quiet=TRUE)
med <- median(lengths)
buff <- med*${len_cutoff}
maximum <- round(med + buff)
cat(maximum)" > ${tmp_dir}/${1}_max.R

gtt-count-bases-per-seq -i ${tmp_dir}/${1}_hits.faa -o ${1}_Num_bps.tmp
cut -f2 ${1}_Num_bps.tmp > ${1}_lengths.tmp
min_len=$(Rscript ${tmp_dir}/${1}_min.R)
max_len=$(Rscript ${tmp_dir}/${1}_max.R)
gtt-filter-seqs-by-length -i ${tmp_dir}/${1}_hits.faa -m $min_len -M $max_len -o ${tmp_dir}/${1}_hits_filtered.tmp > ${1}_filter.out.tmp

cat <(printf "\n    Filtering ${GREEN}${1}${NC} sequences by length...\n") ${1}_filter.out.tmp

rm ${1}_Num_bps.tmp ${1}_lengths.tmp ${1}_filter.out.tmp
