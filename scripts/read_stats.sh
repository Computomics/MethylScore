#!/bin/bash
#$ -S /bin/bash
set -e

sample=$1
filelist=$2
folder=$3
extbin_path=$4
do_dedup=$5
ROI=$6

samtools=$extbin_path/st
if [ ! -e $samtools ]; then >&2 echo "{readstats} Cannot find $samtools"; exit 1; fi
bedtools=$extbin_path/bdt
if [ ! -e $bedtools ]; then >&2 echo "{readstats} Cannot find $bedtools"; exit 1; fi

# split filelist by comma (and newline due to end of command)
# filelist contains only absolute paths to mapped files
files=$(echo $filelist | tr "," "\n")

# calculate total nr of reads in mapping file
totalNR=0
passNR=0
for file in $files; do
    if [ ! -e $file ]; then >&2 echo "{readstats} Cannot find $file"; exit 1; fi
    totalNR=$((totalNR + `$samtools view -c $file`))
    passNR=$((passNR + `$samtools view -c -F 4 -F 256 -F 512 $file`))
done

bamfile=$folder/$sample.passQC.dedup.bam
if [ $do_dedup -eq 0 ]; then
    bamfile=$folder/$sample.passQC.bam
fi
if [ ! -e $bamfile ]; then >&2 echo "{readstats} Cannot find $bamfile"; exit 1; fi



dedupNR=`$samtools view -c -F 256 $bamfile` # passQC already removed unmapped (skip -F 4 here)
uniqNR=`$samtools view -c -F 256 -q 30 $bamfile`

passPERC=`awk -v a=$passNR -v all=$totalNR 'BEGIN{print sprintf("%.1f", -1*(100-100*a/all)) }'`
dedupPERC=`awk -v a=$dedupNR -v all=$totalNR -v u=$passPERC 'BEGIN{print sprintf("%.1f", -1*(100-100*a/all)-u) }'`
uniqPERC=`awk -v a=$uniqNR -v all=$totalNR -v u=$passPERC -v d=$dedupPERC 'BEGIN{print sprintf("%.1f", -1*(100-100*a/all)-u-d) }'`


if [ -e $folder/$sample.read_stats.tsv ]; then
  mv $folder/$sample.read_stats.tsv $folder/$sample.read_stats.tsv.$RANDOM
fi

if [[ $ROI != "" ]]; then
    roiNR=`$bedtools intersect -u -a $bamfile -b $ROI | $samtools view -c -F 256 -q 30 -`
    roiPERC=`awk -v a=$roiNR -v all=$totalNR -v u=$passPERC -v d=$dedupPERC -v r=$uniqPERC 'BEGIN{print sprintf("%.1f", -1*(100-100*a/all)-u-d-r) }'`

    echo -e "#sample\ttotal\tumapped\tdupl\tmultpl\toff-trg\t#ROI\tROI/flt\tROI/tot" > $folder/$sample.read_stats.tsv
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
    ) | paste - - - - - - - - - >> $folder/$sample.read_stats.tsv
else
    echo -e "#sample\ttotal\tumapped\tdupl\tmultpl" > $folder/$sample.read_stats.tsv
    (
    echo $sample
    echo $totalNR
    echo $passPERC
    echo $dedupPERC
    echo $uniqPERC
    ) | paste - - - - - >> $folder/$sample.read_stats.tsv
fi
