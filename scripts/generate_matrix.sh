#!/bin/bash
#$ -S /bin/bash
set -e

stage=$1
bin_path=$2
samplefile=$3
outdir=$4
chr=$5

if [ ! -e $samplefile ]; then >&2 echo "{$stage} Cannot find $samplefile"; exit 1; fi
if [ ! -e "$bin_path/generate_genome_matrix" ]; then
    >&2 echo "{$stage} Cannot find $bin_path/generate_genome_matrix"; exit 1;
fi

echo "### [$stage] Generating genome matrix for chromosome $chr"

$bin_path/generate_genome_matrix -s $samplefile -i mxX -o $outdir/genome_matrix.$chr.tsv -t $outdir

touch $outdir/done
