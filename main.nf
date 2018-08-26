#!/usr/bin/env nextflow
/*
========================================================================================
                                M e t h y l S c o r e
========================================================================================
 Nextflow implementation of Computomics' MethylScore Pipeline
 #### Homepage / Documentation
 https://github.com/Gregor-Mendel-Institute/MethylScore-nf
 #### Author
 Patrick Hüther <patrick.huether@gmi.oeaw.ac.at>
----------------------------------------------------------------------------------------
*/


/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = "0.1.13.2-nf"

// Configurable variables
params.CLUSTER_PROJECT = "becker_common"
params.GENOME = "/lustre/scratch/datasets/TAIR/9/fasta/TAIR9.fa"
params.IGV = 0
params.CLUSTER_MIN_METH = 20
params.CLUSTER_MIN_METH_DIFF = 20
params.DESERT_SIZE = 100
params.DMR_MIN_C = 10
params.DMR_MIN_COV = 3
params.DO_DEDUP = 1

params.FDR_CUTOFF = 0.05
params.FORCE_RERUN = 0
params.HDMR_FOLD_CHANGE = 3
params.HUMAN = 0
params.IGNORE_FIRST_BP = 3
params.IGNORE_LAST_BP = 1
params.MERGE_DIST = 30
params.MIN_COVERAGE = 1
params.MIN_QUAL = 30
params.MR_BATCH_SIZE = 500
params.MR_FREQ_CHANGE = 20
params.MR_FREQ_DISTANCE = 30
params.MR_MIN_C = 20
params.PROJECT_FOLDER = "./results"
params.REMOVE_INTMED_FILES = 0
params.ROI = ""
params.SLIDING_WINDOW_SIZE = 0
params.SLIDING_WINDOW_STEP = 0
params.STATISTICS = 1
params.TRIM_METHRATE = 10

params.SAMPLE_SHEET = "./samplesheet.tsv"
params.SCRIPT_PATH = "scripts"
params.BIN_PATH = "bin"
params.EXTBIN_PATH = "bin_ext"

params.DEBUG = false

// ToDo check parameters for sanity
//if( params.CLUSTER_MIN_METH ){
//    println(params.CLUSTER_MIN_METH.class)
//} else {
//    println(params.CLUSTER_MIN_METH.class)
//}

log.info "=================================================="
log.info " MethylScore ${version}"
log.info "=================================================="
log.info "Reference genome      : ${params.GENOME}"
log.info "Project               : ${params.CLUSTER_PROJECT}"
log.info "Current home          : $HOME"
log.info "Current user          : $USER"
log.info "Current path          : $PWD"
log.info "Script dir            : $baseDir"
log.info "Working dir           : $workDir"
log.info "Output dir            : ${params.PROJECT_FOLDER}"
log.info "---------------------------------------------------"
log.info "IGV output            : ${params.IGV}"
log.info "---------------------------------------------------"
log.info "BIN_PATH              : ${params.BIN_PATH}"
log.info "EXTBIN_PATH           : ${params.EXTBIN_PATH}"
log.info "SCRIPT_PATH           : ${params.SCRIPT_PATH}"
log.info "---------------------------------------------------"
log.info "CLUSTER_MIN_METH      : ${params.CLUSTER_MIN_METH}"
log.info "CLUSTER_MIN_METH_DIFF : ${params.CLUSTER_MIN_METH_DIFF}"
log.info "CLUSTER_PROJECT       : ${params.CLUSTER_PROJECT}"
log.info "DESERT_SIZE           : ${params.DESERT_SIZE}"
log.info "DMR_MIN_C             : ${params.DMR_MIN_C}"
log.info "DMR_MIN_COV           : ${params.DMR_MIN_COV}"
log.info "DO_DEDUP              : ${params.DO_DEDUP ? "Yes" : "No"}"
log.info "FDR_CUTOFF            : ${params.FDR_CUTOFF}"
log.info "FORCE_RERUN           : ${params.FORCE_RERUN}"
log.info "HDMR_FOLD_CHANGE      : ${params.HDMR_FOLD_CHANGE}"
log.info "HUMAN                 : ${params.HUMAN ? "Yes" : "No"}"
log.info "IGNORE_FIRST_BP       : ${params.IGNORE_FIRST_BP}"
log.info "IGNORE_LAST_BP        : ${params.IGNORE_LAST_BP}"
log.info "MERGE_DIST            : ${params.MERGE_DIST}"
log.info "MIN_COVERAGE          : ${params.MIN_COVERAGE}"
log.info "MIN_QUAL              : ${params.MIN_QUAL}"
log.info "MR_BATCH_SIZE         : ${params.MR_BATCH_SIZE}"
log.info "MR_FREQ_CHANGE        : ${params.MR_FREQ_CHANGE}"
log.info "MR_FREQ_DISTANCE      : ${params.MR_FREQ_DISTANCE}"
log.info "MR_MIN_C              : ${params.MR_MIN_C}"
log.info "PROJECT_FOLDER        : ${params.PROJECT_FOLDER}"
log.info "REMOVE_INTMED_FILES   : ${params.REMOVE_INTMED_FILES ? "Yes" : "No"}"
log.info "ROI                   : ${params.ROI}"
log.info "SAMPLE_SHEET          : ${params.SAMPLE_SHEET}"
log.info "SLIDING_WINDOW_SIZE   : ${params.SLIDING_WINDOW_SIZE}"
log.info "SLIDING_WINDOW_STEP   : ${params.SLIDING_WINDOW_STEP}"
log.info "STATISTICS            : ${params.STATISTICS ? "Yes" : "No"}"
log.info "TRIM_METHRATE         : ${params.TRIM_METHRATE}"
log.info "---------------------------------------------------"
log.info "Config Profile : ${workflow.profile}"
log.info "=================================================="

reference = file(params.GENOME)
samplesheet = file(params.SAMPLE_SHEET)

/*
 * Create a channel for the tsv file containing samples
 */

samples = Channel.from(samplesheet.readLines())
                 .map{line ->
                  list = line.split()
                  sampleID = list[0]
                  seqType = list[1]
                  bamFile = list[2]
                  [ sampleID, seqType, file(bamFile) ]
                  }

/*
 * Create a channel for the reference genome and split it by chromosome
 */

refSplit = Channel.fromPath(reference)
	       .splitFasta( record: [id: true, text: true] )

process MethylScore_filterQC {
    tag "$bamFile"
    publishDir "${params.PROJECT_FOLDER}/01mappings", mode: 'copy'

    input:
    set val(sampleID), val(seqType), file(bamFile) from samples

    output:
    set val(sampleID), val(seqType), file('*passQC.bam') into passQC
 
    script:
    """
    samtools view -bh -F 0x200 -F 0x4 -o ${bamFile.baseName}.passQC.bam ${bamFile}
    """
}

process MethylScore_mergeReplicates {
    tag "$sampleID: ${bamFile.collect().size()} replicate(s)"
    publishDir "${params.PROJECT_FOLDER}/01mappings/${sampleID}", mode: 'copy'

    input:
    set val(sampleID), val(seqType), file(bamFile) from passQC.groupTuple()

    output:
    set val(sampleID), val(seqType), file('*passQC.bam') into merged

    script:
    if( bamFile.toList().size() != 1 )
       """
       mkdir tmp

       java -Xmx1024m -Xms256m -Djava.io.tmpdir=tmp -XX:ParallelGCThreads=1 -jar ${PICARD_PATH}/picard.jar \\
        MergeSamFiles \\
    	  I=${bamFile.join(' I=')} \\
    	  O=${sampleID}.passQC.bam \\
    	  USE_THREADING=false
       """
    else
       """
       mv ${bamFile} ${sampleID}.passQC.bam
       """
}

if(params.DO_DEDUP) {

 process MethylScore_deduplicate {
    tag "$sampleID"
    publishDir "${params.PROJECT_FOLDER}/01mappings/${sampleID}", mode: 'copy'

    input:
    set val(sampleID), val(seqType), file(bamFile) from merged

    output:
    set val(sampleID), val(seqType), file('*passQC.dedup.bam') into (dedup, stats)
    file ('dedup.metrics.txt')

    script:
    """
    mkdir tmp
    java -Xmx1024m -Xms256m -Djava.io.tmpdir=tmp -XX:ParallelGCThreads=1 -jar ${PICARD_PATH}/picard.jar \\
      MarkDuplicates \\
        I=${sampleID}.passQC.bam \\
        O=${sampleID}.passQC.dedup.bam \\
        METRICS_FILE=dedup.metrics.txt
        REMOVE_DUPLICATES=true \\
        MAX_FILE_HANDLES=1 \\
        TMP_DIR=\$PWD \\
	VALIDATION_STRINGENCY=LENIENT
    """
 } 
} else {

dedup = merged

}

process MethylScore_readStatistics {
    tag "$sampleID"
    publishDir "${params.PROJECT_FOLDER}/01mappings/${sampleID}", mode: 'copy'

    input:
    set val(sampleID), val(seqType), file(bamFile) from stats

    output:
    file ('*') into readStats

    when:
    params.STATISTICS == 1
   
    script:
    """
    ${baseDir}/${params.SCRIPT_PATH}/read_stats.sh \\
		${sampleID} \\
                ${bamFile} \\
                ${params.ROI}

    if [[ ! -z "${params.ROI}" ]]; then
       ${baseDir}/${params.SCRIPT_PATH}/cov_stats.sh \\
                ${sampleID} \\
                ${bamFile} \\
                ${params.ROI}

    fi
    """
}

process MethylScore_splitBams {
    tag "$sampleID:$chromosome.id"
    publishDir "${params.PROJECT_FOLDER}/01mappings/${sampleID}/split/${chromosome.id}", mode: 'copy'

    input:
    set val(sampleID), val(seqType), file(bamFile) from dedup
    each chromosome from refSplit

    output:
    set val(sampleID), val(seqType), file('*bam'), val("${chromosome.id}"), val("${chromosome.text}") into chrSplit

    script:
    """
    samtools index ${bamFile}
    samtools view -b ${bamFile} ${chromosome.id} > ${bamFile.baseName}_${chromosome.id}.bam
    """
}

process MethylScore_callConsensus {
    tag "$sampleID:$chromosome"
    publishDir "${params.PROJECT_FOLDER}/02consensus/${sampleID}/${chromosome}", mode: 'copy'

    input:
    set val(sampleID), val(seqType), file(splitBam), val(chromosome), file('reference.fa') from chrSplit
//    each context from (['CG','CHG','CHH']) // TODO: check whether parallelizing could be beneficial for performance here

    output:
    set val(sampleID), file('*allC.output'), val(chromosome) into consensus

    script:
    flagW = ( seqType[0] == "PE" ? "flagW=99,147 flagC=83,163" : "flagW=0 flagC=16" )
    """
    MethylExtract.pl \\
     seq='.' \\
     inDir='.' \\
     peOverlap=Y \\
     delDup=N \\
     minQ=${params.MIN_QUAL} \\
     methNonCpGs=0 \\
     varFraction=0.01 \\
     maxPval=0.01 \\
     p=1 \\
     chromDiv=100 \\
     memNumReads=1000 \\
     FirstIgnor=${params.IGNORE_FIRST_BP} \\
     LastIgnor=${params.IGNORE_LAST_BP} \\
     minDepthMeth=1 \\
     context=ALL \\
     bedOut=N wigOut=N \\
     ${flagW}

    for context in CG CHG CHH; do
     sort -k1,1 -k2,2n \$context.output |
     awk -vi=\$context -vs=$sampleID '\$0!~/^#/{OFS="\\t"; \$1=s "\\t" \$1; \$3=i "." \$3; print \$0}' > ${sampleID}.${chromosome}.\$context.output.tmp;
    done
    
    sort -m -k2,2 -k3,3n *.output.tmp > ${sampleID}.${chromosome}.allC.output
    """
}

process MethylScore_chromosomalmatrix {
    tag "$chromosome"
    publishDir "${params.PROJECT_FOLDER}/03matrix", mode: 'copy'

    input:
    set val(sampleID), file(allC), val(chromosome) from consensus.groupTuple(by: 2, sort: 'true')

    output:
    set val(chromosome), file("${chromosome}.genome_matrix.tsv") into splitMatrix
    val(sampleID) into idx

    script:
    """
    idx=4 
    for i in ${sampleID.join(' ')}; do
    	echo -e \$i'\t'\$i'.'${chromosome}'.allC.output' >> samples.txt
    done

    generate_genome_matrix -s samples.txt -i mxX -o ${chromosome}.genome_matrix.tsv
    """
}

splitMatrix
 .transpose()
 .toSortedList({a, b -> a[0] <=> b[0]})
 .flatten()
 .buffer(size:1, skip:1)
 .set { matrix }

process MethylScore_genomematrix {
    tag "$matrixWG"
    publishDir "${params.PROJECT_FOLDER}/03matrix", mode: 'copy'

    input:
//    file(matrixWG) from splitMatrix.collectFile(name: 'genome_matrix.tsv')
    file(matrixWG) from matrix.collect()

    output:
    file('genome_matrix.tsv') into (matrixWG_MRs, matrixWG_igv, matrixWG_DMRs)

    script:
    """
    cat ${matrixWG} | sed '1!{/^#/d;}' > genome_matrix.tsv
    """
//    sed -i '1!{/^#/d;}' ${matrixWG} > genome_matrix.tsv
}

idx
 .flatten()
 .unique()
 .toList()
 .into { indexedSamples_MRs; indexedSamples_splitting }

process MethylScore_callMRs {
    tag "$sampleID"
    publishDir "${params.PROJECT_FOLDER}/04MRs/${sampleID[0]}", mode: 'copy'

    input:
    file(matrixWG) from matrixWG_MRs
    each sampleID from indexedSamples_MRs.withIndex()

    output:
    file('*MRs.bed') into (MRs_igv, MRs_splitting)
    file('*.params') into hmmparams
    file('*.tsv') into MRstats

    script:
    HUMAN = ( params.HUMAN != 0 ? "-human" : "" )
    MIN_C = ( params.MR_MIN_C != 0 ? "-n ${params.MR_MIN_C}" : "" )
    """
    hmm_mrs \\
     -x ${sampleID[1]+1} \\
     -y ${sampleID[0]} \\
     -c ${params.MIN_COVERAGE} \\
     -o ${sampleID[0]}.MRs.bed \\
     -d ${params.DESERT_SIZE} \\
     -i 30 \\
     -m ${params.MERGE_DIST} \\
     -t ${params.TRIM_METHRATE/100} \\
     -p hmm.params \\
     ${HUMAN} \\
     ${MIN_C} \\
     ${matrixWG}

     echo -e "${sampleID[0]}\t" \\
        \$(cat ${sampleID[0]}.MRs.bed | wc -l)"\t" \\
        \$(awk -v OFS=\"\t\" '{sum+=\$3-\$2+1}END{print sum, sprintf("%.0f", sum/NR)}' ${sampleID[0]}.MRs.bed) \\
        > MR_stats.tsv
    """
}

process MethylScore_igv {
    tag "batchsize:${params.MR_BATCH_SIZE}"
    publishDir "${params.PROJECT_FOLDER}/igv", mode: 'copy'

    input:
    file(bed) from MRs_igv.collect()
    file(matrixWG) from matrixWG_igv

    output:
    file('*') into igv

    when:
    params.IGV == true
 
    script:
    """
    sort -m -k1,1 -k2,2g -k3,3g ${bed} > MRs.merged.bed
    python ${baseDir}/${params.SCRIPT_PATH}/matrix2igv.py -i ${matrixWG} -m MRs.merged.bed -o methinfo.igv
    """
}

process MethylScore_splitMRs {
    tag "batchsize:${params.MR_BATCH_SIZE}"
    publishDir "${params.PROJECT_FOLDER}/05DMRs/batches", mode: 'copy'

    input:
    val(sampleID) from indexedSamples_splitting.collect()
    file(MRfile) from MRs_splitting.collect()

    output:
    file('MRbatch*') into MRchunks mode flatten
    file('samples.tsv') into (samples_callDMRs, samples_mergeDMRs)

    script:
    """
    idx=4
    for i in ${sampleID.join(' ')}; do
    	echo -e \$i'\t'\$idx'\t'\$i'.MRs.bed' >> samples.tsv
        (( idx ++ ))
    done

    split_MRfile samples.tsv MRbatch ${params.MR_BATCH_SIZE}
    """
}

process MethylScore_callDMRs {
    tag "$chunk"
    publishDir "${params.PROJECT_FOLDER}/05DMRs", mode: 'copy'

    input:
    file(matrixWG) from matrixWG_DMRs
    file(samples) from samples_callDMRs
    each file(chunk) from MRchunks

    output:
    file('*/*.bed') optional true into bedFiles
    file('*/*.dif') optional true into segmentFiles

    script:
    """
    mkdir ${chunk}.out

    dmrs \\
     -s ${samples} \\
     -r ${chunk} \\
     -m ${matrixWG} \\
     -p ${params.MR_FREQ_CHANGE} \\
     -i ${params.CLUSTER_MIN_METH_DIFF} \\
     -j ${params.CLUSTER_MIN_METH} \\
     -v ${params.DMR_MIN_COV} \\
     -n ${params.DMR_MIN_C} \\
     -w ${params.SLIDING_WINDOW_SIZE} \\
     -x ${params.SLIDING_WINDOW_STEP} \\
     -z 1 \\
     -B ${baseDir}/${params.BIN_PATH}/betabin_model \\
     -Y ${baseDir}/${params.SCRIPT_PATH}/pv2qv.py \\
     --no-post-process \\
     -o ${chunk}.out
    """
}

process MethylScore_mergeDMRs {
    tag "$segments"
    publishDir "${params.PROJECT_FOLDER}/05DMRs", mode: 'copy'

    input:
    file(segments) from segmentFiles.collectFile(name:'segments.dif')
    file(samples) from samples_mergeDMRs

    output:
    file('*') into DMRs

    script:
    """
    merge_DMRs \\
     ${samples} \\
     ${segments} \\
     . \\
     python \\
     ${baseDir}/${params.SCRIPT_PATH}/pv2qv.py \\
     ${params.FDR_CUTOFF} \\
     ${params.CLUSTER_MIN_METH} \\
     ${params.DMR_MIN_C} \\
     ${params.HDMR_FOLD_CHANGE} \\
     ${params.REMOVE_INTMED_FILES}
    """
}
