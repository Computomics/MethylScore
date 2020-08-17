#!/bin/bash
#$ -S /bin/bash
set -e

stage=$1
script_path=$2
matrix=$3
outdir=$4
mrs=$5

if [ ! -e "$script_path/matrix2igv.py" ]; then >&2 echo "{igv} Cannot find $script_path/matrix2igv.py"; exit 1; fi
if [ ! -e "$matrix" ]; then >&2 echo "{igv} Cannot find $matrix"; exit 1; fi
if [[ "$mrs" == "" ]]; then >&2 echo "{igv} Cannot find any MRs"; exit 1; fi

echo -e "### [$stage] Running igv\n"

if [ ! -e $outdir ]; then mkdir -p $outdir; fi
sort -m -k1,1 -k2,2g -k3,3g $mrs > $outdir/MRs.merged.bed

python $script_path/matrix2igv.py -i $matrix -m $outdir/MRs.merged.bed -o $outdir/methinfo.igv

# TODO convert to tdf?
# igvtools toTDF --fileType igv -t . " + args.output + " " + args.output.replace(".igv", ".tdf")

touch $outdir/done
