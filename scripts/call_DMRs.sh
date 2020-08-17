#!/bin/bash
#$ -S /bin/bash
set -e

stage=$1
command=$2
outdir=$3

echo -e "### [$stage] Calling DMRs\n"

$command

touch $outdir/done
