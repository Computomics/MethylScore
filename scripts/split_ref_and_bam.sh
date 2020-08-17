#!/bin/bash
#$ -S /bin/bash
set -e

stage=$1
bamfile=$2
reffile=$3
rerun=$4
bamtools_split=$5

echo "### [$stage] Splitting ref and bam"

# disentagnle bamtools command (replace ',' by ' '):
bamtools_cmd=$(echo $bamtools_split | tr "," " ")
echo -e "### [split_ref_and_bam] Sanity check: bamtools: '$bamtools_cmd'"

if [ ! -e $bamfile ]; then >&2 echo "{splitbam} Cannot find $bamfile"; exit 1; fi

dir=$(dirname $bamfile)
file=$(basename $bamfile)
if [ ! -e "$dir/split" ]; then mkdir "$dir/split"; fi

### STEP 1: split bamfile

    echo "### [$stage] Splitting bamfile by chromosome"

    $bamtools_cmd -reference -in $bamfile

    for f in $dir/${file%.bam}.REF*; do
        chr=${f#*REF_}
        chr=${chr%.bam}
        if [ ! -e $dir/split/$chr ]; then mkdir $dir/split/$chr; fi
        mv $dir/${file%.bam}.REF_$chr.bam $dir/split/$chr/
    done


### STEP 2: split reference

    if [ ! -e $reffile ]; then >&2 echo "{splitbam} Cannot find $reffile"; exit 1; fi
    echo "### [$stage] Splitting reference by chromosome"

    while read line
    do
        if [[ ${line:0:1} == '>' ]]
        then
            chr=${line%% *}
            if [ ! -e $dir"/split/"${chr#>} ]; then
                outfile="/dev/null"
            else
                outfile=$dir"/split/"${chr#>}/${chr#>}".fa"
            fi

            echo $line > $outfile
        else
            echo $line >> $outfile
        fi
    done < $reffile
