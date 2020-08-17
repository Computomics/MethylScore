#!/bin/bash
#$ -S /bin/bash
set -e

stage=$1
bin_path=$2
sample=$3
sample_idx=$4
cov=$5
desertsize=$6
mergedist=$7
trimlevel=$8
matrix=$9
outdir=${10}
human=${11}
min_c=${12}
params=${13}

if [ ! -e "$bin_path/hmm_mrs" ]; then >&2 echo "{$stage} Cannot find $bin_path/hmm_mrs"; exit 1; fi
if [ ! -e "$matrix" ]; then >&2 echo "{$stage} Cannot find $matrix"; exit 1; fi

echo -e "### [$stage] Calling MRs\n"

if [ ! -e $outdir ]; then mkdir -p $outdir; fi

min_c_param=""
if [ $min_c -gt 0 ]; then min_c_param="-n $min_c"; fi

human_param=""
if [ $human -eq 1 ]; then human_param="-human"; fi

params_param=""
if [ $params != "" ]; then params_param="-P $params"; fi

# call MRs:
$bin_path/hmm_mrs \
    -x $sample_idx \
    -y $sample \
    -c $cov \
    -o $outdir/MRs.bed \
    -d $desertsize \
    -i 30 \
    -m $mergedist \
    -t $trimlevel \
    -p $outdir/hmm.params \
    $human_param \
    $min_c_param \
    $params_param \
    $matrix

# calculate little statistics:
if [ -e $outdir/MR_stats.tsv ]; then
  mv $outdir/MR_stats.tsv \
     $outdir/MR_stats.tsv.$RANDOM
fi

sample=${outdir##*/}

echo -e "$sample\t" \
        `cat $outdir/MRs.bed|wc -l`"\t" \
        `awk -v OFS=\"\t\" '{sum+=$3-$2+1}END{print sum, sprintf("%.0f", sum/NR)}' $outdir/MRs.bed` \
        > $outdir/MR_stats.tsv

touch $outdir/done
