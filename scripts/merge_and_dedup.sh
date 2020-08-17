#!/bin/bash
#$ -S /bin/bash
set -e

stage=$1
mapping_path=$2
sample=$3
samplelist=$4
reference=$5
extbin_path=$6
rerun=$7
rm_intmd=$8
max_mem=$9"m"
bamtools_split=${10}
samtools_cmd=${11}
stats=${12}
bin_path=${13}
do_dedup=${14}
rois=${15}


if [ ! -e "$extbin_path/picard.jar" ]; then >&2 echo "{dedup} Cannot find $extbin_path/picard.jar"; exit 1; fi

echo "### [$stage] Merge and dedup"

outpath=$mapping_path/$sample
if [ ! -e "$outpath" ]; then mkdir -p $outpath; fi



### STEP --1-- QC filter for each replicate

# split samplelist by comma (and newline due to end of command)
# samplelist contains only absolute paths to mapped files
samples=$(echo $samplelist | tr "," "\n")

# Iterate over potential replicates:
file_str=""
file_list=""
nr=0
for infile in $samples; do

    if [ ! -e $infile ]; then >&2 echo "{dedup} Cannot find $infile"; exit 1; fi

    filename=$(basename $infile)
    outfile=$outpath"/"${filename/.bam/.passQC.bam}

    ### Filter out unmapped reads and reads having the flag "not passing platform/vendor QC"
    if [ ! -e "$outpath/qc.done" -o $rerun -eq 1 ]; then
        echo "### [$stage] PassQC"
        $samtools_cmd view -bh -F 0x200 -F 0x4 -o $outfile $infile
        #if [ $rm_intmd -eq 1 ]; then rm $infile; fi
    fi

    # generate a string: "I=sample1 I=sample2 ..."
    file_str=$file_str" I="$outfile
    file_list=$file_list" "$outfile
	  nr=$((nr+1))

done




### STEP --2-- Merge technical replicates if applicable

if [ $nr -gt 1 ]; then

    # if there is more than one technical replicate, merge bam files:
    echo "### [$stage] Merge replicates"

	  java -Xmx$max_mem -Xms$max_mem -XX:ParallelGCThreads=1 -jar $extbin_path/picard.jar \
	    MergeSamFiles \
	      $file_str \
	      O=$outpath/${sample}.passQC.bam \
	      USE_THREADING=false

else

    filename=$(basename $samplelist)
    if [ ${filename%.bam} != $sample ]; then mv $outpath/${filename/.bam/.passQC.bam} $outpath/$sample.passQC.bam; fi

fi



### STEP --3-- Remove duplicate reads

if [ $do_dedup -eq 1 ]; then
    echo "### [$stage] De-duplicating"

    java -Xmx$max_mem -Xms$max_mem -XX:ParallelGCThreads=1 -jar $extbin_path/picard.jar \
      MarkDuplicates \
        INPUT=$outpath/${sample}.passQC.bam \
        OUTPUT=$outpath/${sample}.passQC.dedup.bam \
        METRICS_FILE=$outpath/dedup.metrics.txt \
        REMOVE_DUPLICATES=true \
        MAX_FILE_HANDLES=1 \
        TMP_DIR=$outpath \
        VALIDATION_STRINGENCY=LENIENT

    bamfile=$outpath/${sample}.passQC.dedup.bam
else
    bamfile=$outpath/${sample}.passQC.bam
fi

    echo "### [$stage] Indexing bamfile"
    $samtools_cmd index $bamfile


### do statistcs:
if [ $stats -eq 1 ]; then
    if [ ! -e "$bin_path/read_stats.sh" ]; then >&2 echo "{dedup} Cannot find $bin_path/read_stats.sh"; exit 1; fi
    echo "### [$stage] Collecting read statistics"
    bash $bin_path/read_stats.sh \
        $sample \
        $samplelist \
        $outpath \
        $extbin_path \
        $do_dedup \
        $rois

    if [[ $rois != "" ]]; then
        if [ ! -e "$bin_path/cov_stats.sh" ]; then >&2 echo "{dedup} Cannot find $bin_path/cov_stats.sh"; exit 1; fi
        echo "### [$stage] Collecting coverage statistics"
        bash $bin_path/cov_stats.sh \
            $sample \
            $outpath \
            $extbin_path \
            $do_dedup \
            $rois
    fi
fi

### possibly remove intermediate files:
if [ $rm_intmd -eq 1 -a $do_dedup -eq 1 ]; then
  for file in $outpath/*passQC.bam; do
    rm $file
  done
fi



### STEP --4-- Split bam and ref file for better parallelization

if [ ! -e "$bin_path/split_ref_and_bam.sh" ]; then >&2 echo "{dedup} Cannot find $bin_path/split_ref_and_bam.sh"; exit 1; fi
bash $bin_path/split_ref_and_bam.sh \
    $stage \
    $bamfile \
    $reference \
    $rerun \
    $bamtools_split


touch $outpath/done
