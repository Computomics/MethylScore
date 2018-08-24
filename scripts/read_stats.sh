#!/bin/bash
set -e

sample=$1
bamfile=$2
ROI=$3

if [ ! -e $bamfile ]; then >&2 echo "{readstats} Cannot find $bamfile"; exit 1; fi

# calculate total nr of reads in mapping file
totalNR=$((totalNR + `samtools view -c $bamfile`))
passNR=$((passNR + `samtools view -c -F 4 -F 256 -F 512 $bamfile`))
dedupNR=`samtools view -c -F 256 $bamfile` # passQC already removed unmapped (skip -F 4 here)
uniqNR=`samtools view -c -F 256 -q 30 $bamfile`
passPERC=`awk -v a=$passNR -v all=$totalNR 'BEGIN{print sprintf("%.1f", -1*(100-100*a/all)) }'`
dedupPERC=`awk -v a=$dedupNR -v all=$totalNR -v u=$passPERC 'BEGIN{print sprintf("%.1f", -1*(100-100*a/all)-u) }'`
uniqPERC=`awk -v a=$uniqNR -v all=$totalNR -v u=$passPERC -v d=$dedupPERC 'BEGIN{print sprintf("%.1f", -1*(100-100*a/all)-u-d) }'`

if [ -e $sample.read_stats.tsv ]; then
  mv $sample.read_stats.tsv $folder/$sample.read_stats.tsv.$RANDOM
fi

if [[ $ROI != "" ]]; then
    roiNR=`bedtools intersect -u -a $bamfile -b $ROI | samtools view -c -F 256 -q 30 -`
    roiPERC=`awk -v a=$roiNR -v all=$totalNR -v u=$passPERC -v d=$dedupPERC -v r=$uniqPERC 'BEGIN{print sprintf("%.1f", -1*(100-100*a/all)-u-d-r) }'`

    echo -e "#sample\ttotal\tumapped\tdupl\tmultpl\toff-trg\t#ROI\tROI/flt\tROI/tot" > $sample.read_stats.tsv
    (
    echo $sample
    echo $totalNR
    echo $passPERC
    echo $dedupPERC
    echo $uniqPERC
    echo $roiPERC
    echo $roiNR
    awk -v a=$roiNR -v b=$uniqNR 'BEGIN{print sprintf("%.1f", 100*a/b) }'
    awk -v a=$roiNR -v b=$totalNR 'BEGIN{print sprintf("%.1f", 100*a/b) }'
    ) | paste - - - - - - - - - >> $sample.read_stats.tsv
else
    echo -e "#sample\ttotal\tumapped\tdupl\tmultpl" > $sample.read_stats.tsv
    (
    echo $sample
    echo $totalNR
    echo $passPERC
    echo $dedupPERC
    echo $uniqPERC
    ) | paste - - - - - >> $sample.read_stats.tsv
fi
