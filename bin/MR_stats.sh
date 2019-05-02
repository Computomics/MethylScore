#!/usr/bin/env bash
set -e

sample=$1
bedfile=$2

echo -e "#sampleID\tMRcount\tgenomespace\tavg_length" > $sample.MR_stats.tsv
echo -e "$sample\t$(cat $bedfile | wc -l)\t$(awk -v OFS='\t' '{sum+=$3-$2+1}END{print sum, sprintf("%.0f", sum/NR)}' $bedfile)" >> $sample.MR_stats.tsv
