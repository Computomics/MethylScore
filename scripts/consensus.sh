#!/bin/bash
#$ -S /bin/bash
set -e

stage=$1
extbin_path=$2
sample=$3
libtype=$4
minqual=$5
ignore_first=$6
ignore_last=$7
indir=$8
outdir=$9
rm_intmd=${10}
samtools_cmd=${11}


if [ ! -e "$extbin_path/MethylExtract.pl" ]; then >&2
    echo "[$stage] Cannot find $extbin_path/MethylExtract.pl"
    exit 1
fi

if [ ! -e $indir/tmp ]; then mkdir $indir/tmp; fi

echo -e "### [$stage] Running consensus\n"

sam_flags=""
if [ $libtype == "SE" ]; then
    sam_flags="flagW=0 flagC=16"
elif [ $libtype == "PE" ]; then
    sam_flags="flagW=99,147 flagC=83,163"
fi


### STEP 1: run consensus

perl $extbin_path/MethylExtract.pl \
    seq=$indir/ \
    inDir=$indir/ \
    peOverlap=Y \
    delDup=N \
    minQ=$minqual \
    methNonCpGs=0 \
    varFraction=0.01 \
    maxPval=0.01 \
    p=1 \
    chromDiv=100 \
    memNumReads=1000 \
    FirstIgnor=$ignore_first \
    LastIgnor=$ignore_last \
    minDepthMeth=1 \
    context=ALL \
    bedOut=N wigOut=N \
    $sam_flags \
    outDir=$outdir \
    samtools=$samtools_cmd


### STEP 2: sort and merge output data

#sort -k1,1 -k2,2n $outdir/*.output | grep -v "^#" > $outdir/allC.output

for i in CG CHG CHH; do
    if [ ! -e "$outdir/$i.output" ]; then >&2 echo "[$stage] Cannot find $outdir/$i.output"; exit 1; fi
    sort -k1,1 -k2,2n -T $indir/tmp $outdir/$i.output |
        awk -vi=$i -vs=$sample '$0!~/^#/{OFS="\t"; $1=s "\t" $1; $3=i "." $3; print $0}' > $outdir/$i.output.tmp
done

sort -m -k2,2 -k3,3n -T $indir/tmp $outdir/*.output.tmp > $outdir/allC.output
rm $outdir/*.output.tmp
if [ $rm_intmd -eq 1 ]; then rm -rf $outdir/C*output $indir/*; fi

# make sure the output file exists and is copied completely:
while [ ! -e $outdir/allC.output ]; do
  sleep 2;
done
filesize=$(stat -c%s "$outdir/allC.output")
sleep 1
while [ $(stat -c%s "$outdir/allC.output") -ne $filesize ]; do
  filesize=$(stat -c%s "$outdir/allC.output")
  sleep 2;
done

rm -rf $indir/tmp
touch $outdir/done
