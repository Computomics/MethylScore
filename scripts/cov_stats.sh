#!/bin/bash
set -e

sample=$1
bamfile=$2
rois=$3

if [ ! -e $bamfile ]; then >&2 echo "Cannot find $bamfile"; exit 1; fi
if [ ! -e $rois ]; then >&2 echo "Cannot find $rois"; exit 1; fi

bp_ROIs=`awk '{sum += $3-$2} END{print sum}' $rois`

bedtools genomecov -bg -ibam <( samtools view -bh -q 30 -F 4 -F 256 $bamfile ) | \
tee sample.cov.bed | \
awk -v s=$sample '{sum+=($3-$2)*$4; pos+=$3-$2} END{print s "\t" sprintf("%.1f", sum/pos)}' \
> $sample.cov.avg


if [ -e $sample.cov_stats.tsv ]; then
  mv $sample.cov_stats.tsv $sample.cov_stats.tsv.$RANDOM
fi


echo -e "#sample\tON#pos\tON%pos\tONavgcv\tOFF#pos\tOFFavgc\tenrich" > $sample.cov_stats.tsv
#### Output format columns ####
#                 ON/OFF-TARGET
#sample name
#covered_bp       ON
#covered_%        ON
#avg_cov_mincov1  ON
#covered_bp       OFF
#avg_cov_mincov1  OFF
#enrichment factor (ON/ON+OFF)
(
echo $sample
bedtools intersect -wo -a $rois -b $sample.cov.bed | awk -v bp=$bp_ROIs '{ if($2<=$6) start=$6;  if($2>$6) start=$2; if($3<=$7) end=$3; if($3>$7) end=$7; len=end-start; sum+=len*$8; nr_pos+=len; } END{ OFS="\t"; print nr_pos, sprintf("%.0f", 100*nr_pos/bp), sprintf("%.1f", sum/nr_pos) }'
bedtools subtract -a $sample.cov.bed -b $rois | awk '{ len=$3-$2; sum+=len*$4; nr_pos+=len } END{ OFS="\t"; print nr_pos, sprintf("%.1f", sum/nr_pos) }'
) | paste - - - - - - | awk '{ OFS="\t"; if ($6>0) print $1,$2,$3,$4,$5,$6,sprintf("%.1f", $4/$6); else print $0 }' >> $sample.cov_stats.tsv

rm sample.cov.bed
