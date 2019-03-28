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
version = "0.1.14-nf"

// General parameters
params.CLUSTER_PROJECT = "becker_common"
params.GENOME = "/lustre/scratch/datasets/TAIR/9/fasta/TAIR9.fa"
params.SAMPLE_SHEET = false
params.BEDGRAPH = false
params.PROJECT_FOLDER = "./results"

params.HUMAN = false
params.ROI = false

params.SCRIPT_PATH = "scripts"
params.BIN_PATH = "bin"
params.EXTBIN_PATH = "bin_ext"

params.STATISTICS = true
params.IGV = false
params.DO_DEDUP = true

// DMR parameters
params.MR_FREQ_CHANGE = 20
params.MR_FREQ_DISTANCE = 30
params.CLUSTER_MIN_METH = 20
params.CLUSTER_MIN_METH_DIFF = 20
params.SLIDING_WINDOW_SIZE = 0
params.SLIDING_WINDOW_STEP = 0
params.DMR_MIN_C = 10
params.DMR_MIN_COV = 3
params.MR_BATCH_SIZE = 500
params.HDMR_FOLD_CHANGE = 3
params.FDR_CUTOFF = 0.05

// Consensus parameters
params.MIN_QUAL = 30
params.IGNORE_FIRST_BP = 3
params.IGNORE_LAST_BP = 1

// MR parameters
params.MIN_COVERAGE = 1
params.MR_MIN_C = 20
params.DESERT_SIZE = 100
params.TRIM_METHRATE = 10
params.MERGE_DIST = 30
params.MR_PARAMS = false

params.DEBUG = false

// Parameter checks
assert params.SAMPLE_SHEET, "samplesheet.tsv has to be specified!"
assert params.GENOME, "reference genome in fasta format has to be specified!"

assert params.HUMAN instanceof Boolean, "HUMAN must be set to either false (off) or true (on)"
assert params.IGV instanceof Boolean, "IGV must be set to either false (off) or true (on)"
assert params.STATISTICS instanceof Boolean, "STATISTICS must be set to either false (off) or true (on)"
assert params.DO_DEDUP instanceof Boolean, "DO_DEDUP must be set to either false (off) or true (on)"

assert params.MR_FREQ_CHANGE in 0..100, "MR_FREQ_CHANGE must be between 0 and 100!"
assert params.CLUSTER_MIN_METH_DIFF in 0..100, "CLUSTER_MIN_METH_DIFF must be between 0 and 100!"
assert params.CLUSTER_MIN_METH in 0..100, "CLUSTER_MIN_METH must be between 0 and 100!"
assert params.MR_FREQ_DISTANCE instanceof Integer && params.MR_FREQ_DISTANCE >= 0, "MR_FREQ_DISTANCE must be a non-negative integer!"
assert params.SLIDING_WINDOW_SIZE instanceof Integer && params.SLIDING_WINDOW_SIZE >= 0, "SLIDING_WINDOW_SIZE must be a non-negative integer!"
assert params.SLIDING_WINDOW_STEP instanceof Integer && params.SLIDING_WINDOW_STEP >= 0, "SLIDING_WINDOW_STEP must be a non-negative integer!"
assert params.DMR_MIN_C instanceof Integer && params.DMR_MIN_C >= 0, "DMR_MIN_C must be a non-negative integer!"
assert params.DMR_MIN_COV instanceof Integer && params.DMR_MIN_COV >= 0, "DMR_MIN_COV must be a non-negative integer!"
assert params.MR_BATCH_SIZE instanceof Integer && params.MR_BATCH_SIZE >= 0, "MR_BATCH_SIZE must be a non-negative integer!"
assert params.HDMR_FOLD_CHANGE >= 0, "HDMR_FOLD_CHANGE must be a non-negative number!"
assert params.FDR_CUTOFF > 0 && params.FDR_CUTOFF < 1, "FDR_CUTOFF must be between 0 and 1!"

assert params.MIN_QUAL instanceof Integer && params.MIN_QUAL in 1..40, "MIN_QUAL must be between 1 and 40!"
assert params.IGNORE_FIRST_BP instanceof Integer && params.IGNORE_FIRST_BP >= 0, "IGNORE_FIRST_BP must be a non-negative integer!"
assert params.IGNORE_LAST_BP instanceof Integer && params.IGNORE_LAST_BP >= 0, "IGNORE_LAST_BP must be a non-negative integer!"

assert params.MIN_COVERAGE instanceof Integer && params.MIN_COVERAGE > 0, "MIN_COVERAGE must be a non-negative integer!"
assert params.DESERT_SIZE instanceof Integer && params.DESERT_SIZE > 0, "DESERT_SIZE must be a non-negative integer!"
assert params.MERGE_DIST instanceof Integer && params.MERGE_DIST > 0, "MERGE_DIST must be a non-negative integer!"
assert params.MR_MIN_C instanceof Integer, "MR_MIN_C must be an integer! Negative value turns on permutation test."
assert params.TRIM_METHRATE in 0..100, "TRIM_METHRATE must be between 0 and 100!"


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
log.info "DO_DEDUP              : ${params.DO_DEDUP}"
log.info "FDR_CUTOFF            : ${params.FDR_CUTOFF}"
log.info "HDMR_FOLD_CHANGE      : ${params.HDMR_FOLD_CHANGE}"
log.info "HUMAN                 : ${params.HUMAN}"
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
log.info "ROI                   : ${params.ROI}"
log.info "MR_PARAMS             : ${params.MR_PARAMS}"
log.info "SAMPLE_SHEET          : ${params.SAMPLE_SHEET}"
log.info "SLIDING_WINDOW_SIZE   : ${params.SLIDING_WINDOW_SIZE}"
log.info "SLIDING_WINDOW_STEP   : ${params.SLIDING_WINDOW_STEP}"
log.info "STATISTICS            : ${params.STATISTICS}"
log.info "TRIM_METHRATE         : ${params.TRIM_METHRATE}"
log.info "---------------------------------------------------"
log.info "Config Profile : ${workflow.profile}"
log.info "=================================================="

roi_file = params.ROI ? Channel.fromPath("${params.ROI}", checkIfExists: true).collect() : file('null')

hmm_params_file = params.MR_PARAMS ? Channel.fromPath("${params.MR_PARAMS}", checkIfExists: true).collect() : file('null')

/*
 * Create a channel for the reference genome and split it by chromosome
 */

Channel
 .fromPath("${params.GENOME}", checkIfExists: true)
 .splitFasta( record: [id: true, text: true] )
 .set { fasta }

(fasta_MethylExtract, fasta_bedGraph) = !params.BEDGRAPH ? [ fasta, Channel.empty() ] : [ Channel.empty(), fasta ]

/*
 * Create a Channel for the tsv file containing samples
 * Store entries in a map first, then subset to arrays
 */

Channel
  .fromPath("${params.SAMPLE_SHEET}", checkIfExists: true)
  .splitText()
  .map{ line ->
        def list = line.split()
        !params.BEDGRAPH ? [ 'sampleID':list[0], 'seqType':list[1], 'filePath':file(list[2]) ] : [ 'sampleID':list[0], 'filePath':file(list[1]) ]
      }
  .set { inputMap }

def sampleIndex = 1
inputMap
  .map { record -> !params.BEDGRAPH ? [ record.sampleID, record.seqType, record.filePath ] : [ record.sampleID, record.filePath ] }
  .tap { samples }
  .groupTuple()
  .map { record -> tuple( record[0], sampleIndex++ ) }
  .set{ indexedSamples }

(samples_bam, samples_bedGraph) = !params.BEDGRAPH ? [ samples, Channel.empty() ] : [ Channel.empty(), samples ]

/*
 * Start pipeline
 */

process MethylScore_filterQC {
    tag "$bamFile"
    publishDir "${params.PROJECT_FOLDER}/01mappings", mode: 'copy'

    input:
    set val(sampleID), val(seqType), file(bamFile) from samples_bam

    output:
    set val(sampleID), val(seqType), file('*passQC.bam') into passQC
 
    when:
    !params.BEDGRAPH

    script:
    """
    if [[ \$(samtools view -H ${bamFile} | grep unsorted) ]]; then
     samtools sort ${bamFile} | samtools view -bh -F 0x200 -F 0x4 -o ${bamFile.baseName}.passQC.bam
    else
     samtools view -bh -F 0x200 -F 0x4 -o ${bamFile.baseName}.passQC.bam ${bamFile}
    fi
    """
}

process MethylScore_mergeReplicates {
    tag "$sampleID: ${bamFile.collect().size()} replicate(s)"
    publishDir "${params.PROJECT_FOLDER}/01mappings/${sampleID}", mode: 'copy'

    input:
    set val(sampleID), val(seqType), file(bamFile) from passQC.groupTuple()

    output:
    set val(sampleID), val(seqType), file('*merged.passQC.bam') into mergedSamples

    when:
    !params.BEDGRAPH

    script:
    if( bamFile.toList().size() != 1 )
       """
       mkdir tmp
       picard -Xmx${task.memory.toMega()}m -Xms${task.memory.toMega() / 4}m -Djava.io.tmpdir=tmp -XX:ParallelGCThreads=1 \\
        MergeSamFiles \\
        I=${bamFile.join(' I=')} \\
        O=${sampleID}.merged.passQC.bam \\
        USE_THREADING=false
       """
    else
       """
       mv ${bamFile} ${sampleID}.merged.passQC.bam
       """
}


if( params.DO_DEDUP ) {

 process MethylScore_deduplicate {
    tag "$sampleID"
    publishDir "${params.PROJECT_FOLDER}/01mappings/${sampleID}", mode: 'copy'

    input:
    set val(sampleID), val(seqType), file(bamFile) from mergedSamples

    output:
    set val(sampleID), val(seqType), file('*passQC.dedup.bam') into (dedup, stats)
    file ('dedup.metrics.txt')

    when:
    !params.BEDGRAPH

    script:
    """
    mkdir tmp
    picard -Xmx${task.memory.toMega()}m -Xms${task.memory.toMega() / 4}m -Djava.io.tmpdir=tmp -XX:ParallelGCThreads=1 \\
      MarkDuplicates \\
        I=${bamFile} \\
        O=${sampleID}.passQC.dedup.bam \\
        METRICS_FILE=dedup.metrics.txt \\
        REMOVE_DUPLICATES=true \\
        MAX_FILE_HANDLES=1 \\
        TMP_DIR=. \\
        VALIDATION_STRINGENCY=LENIENT
    """
}

} else {

(dedup, stats) = mergedSamples.into(2)

}

process MethylScore_readStatistics {
    tag "$sampleID"
    publishDir "${params.PROJECT_FOLDER}/01mappings/${sampleID}", mode: 'copy'

    input:
    set val(sampleID), val(seqType), file(bamFile) from stats
    file(bed) from roi_file

    output:
    file ('*') into readStats

    when:
    params.STATISTICS && !params.BEDGRAPH
   
    script:
    def REGIONS_FILE = bed.name != 'null' ? "${bed}" : ""

    """
    ${baseDir}/${params.SCRIPT_PATH}/read_stats.sh \\
     ${sampleID} \\
     ${bamFile} \\
     ${REGIONS_FILE}

    ${baseDir}/${params.SCRIPT_PATH}/cov_stats.sh \\
     ${sampleID} \\
     ${bamFile} \\
     ${REGIONS_FILE}
    """
}

process MethylScore_splitBams {
    tag "$sampleID:$chromosome.id"
    publishDir "${params.PROJECT_FOLDER}/01mappings/${sampleID}/split/${chromosome.id}", mode: 'copy'

    input:
    set val(sampleID), val(seqType), file(bamFile) from dedup
    each chromosome from fasta_MethylExtract

    output:
    set val(sampleID), val(seqType), file('*bam'), val("${chromosome.id}"), val("${chromosome.text}") into chrSplit

    when:
    !params.BEDGRAPH

    script:
    """
    samtools index ${bamFile}
    samtools view -b ${bamFile} ${chromosome.id} > ${bamFile.baseName}_${chromosome.id}.bam
    """
}

process MethylScore_callConsensus {
    tag "$sampleID:$chromosomeID"
    publishDir "${params.PROJECT_FOLDER}/02consensus/${sampleID}/${chromosomeID}", mode: 'copy'

    input:
    set val(sampleID), val(seqType), file(splitBam), val(chromosomeID), file("${chromosomeID}.fa") from chrSplit
//    each context from (['CG','CHG','CHH']) // TODO: check whether parallelizing could be beneficial for performance here

    output:
    set val(sampleID), file("${sampleID}.allC.output"), val(chromosomeID), val("${chromosomeID}.fa") into consensus

    when:
    !params.BEDGRAPH

    script:
    def flagW = ( seqType[0] == "PE" ? "flagW=99,147 flagC=83,163 peOverlap=Y" : "flagW=0 flagC=16" )

    """
    MethylExtract.pl \\
     seq='.' \\
     inDir='.' \\
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
     bedOut=N \\
     wigOut=N \\
     ${flagW}

    for context in CG CHG CHH; do
     sort -k1,1 -k2,2n \$context.output |
     awk -vi=\$context -vs=$sampleID '\$0!~/^#/{OFS="\\t"; \$1=s "\\t" \$1; \$3=i "." \$3; print \$0}' > ${sampleID}.${chromosomeID}.\$context.output.tmp;
    done
    
    sort -m -k2,2 -k3,3n *.output.tmp > ${sampleID}.allC.output
    """
}

process MethylScore_mergeContexts {
    tag "$sampleID"
    publishDir "${params.PROJECT_FOLDER}/02consensus", mode: 'copy'

    input:
    set val(sampleID), file(bedGraph) from samples_bedGraph.groupTuple()
    each chromosome from fasta_bedGraph

    output:
    set val(sampleID), file("${sampleID}.merged.bedGraph"), val("${chromosome.id}"), val("${chromosome.text}") into mergedContexts

    when:
    params.BEDGRAPH

    script:
    """
    awk '\$1 == "${chromosome.id}"' ${bedGraph} | sort -k1,1d -k2,2g -T . > ${sampleID}.merged.bedGraph
    """
}

indexedSamples
 .combine(consensus.mix(mergedContexts), by: 0)
 .groupTuple(by: [3,4])
 .set {pile}

process MethylScore_chromosomalmatrix {
    tag "${chromosomeID}"
    publishDir "${params.PROJECT_FOLDER}/03matrix", mode: 'copy'

    input:
    set val(sampleID), val(Index), file(consensus), val(chromosomeID), file("${chromosomeID}.fa") from pile

    output:
    file("${chromosomeID}.genome_matrix.tsv") into splitMatrix
    set val(sampleID), val(Index) into (indexedSamples_splitting, indexedSamples_MRs) mode flatten

    script:
    def inputFormat = params.BEDGRAPH ? "-i bismark -r ${chromosomeID}.fa": "-i mxX"

    """
    paste <(printf "${sampleID.sort{ it.getAt( Index.sort()) }.join('\n')}") \\
          <(printf "${consensus.toList().sort{ it.baseName.getAt( Index.sort() ) }.join('\n')}") > ${chromosomeID}_samples.tsv

    generate_genome_matrix \\
     -s ${chromosomeID}_samples.tsv \\
     ${inputFormat} \\
     -o ${chromosomeID}.genome_matrix.tsv
    """
}

splitMatrix
 .collectFile(name: 'genome_matrix.tsv', keepHeader: true, sort: { it.baseName }, storeDir: "${params.PROJECT_FOLDER}/03matrix")
 .into {matrixWG_MRs; matrixWG_igv; matrixWG_DMRs}

process MethylScore_callMRs {
    tag "${sample.getAt(1)}:${sample.getAt(0)}"
    publishDir "${params.PROJECT_FOLDER}/04MRs/${sample.getAt(0)}", mode: 'copy'

    input:
    file(matrixWG) from matrixWG_MRs
    each sample from indexedSamples_MRs.unique()
    file(parameters) from hmm_params_file

    output:
    file("${sample.getAt(0)}.MRs.bed") into (MRs_igv, MRs_splitting)
    file("${sample.getAt(0)}.hmm_params") optional true into hmm_params
    file("${sample.getAt(0)}.MR_stats.tsv") into MRstats

    script:
    def HUMAN = params.HUMAN ? "-human" : ""
    def MIN_C = params.MR_MIN_C > 0 ? "-n ${params.MR_MIN_C}" : "-n -1"
    def HMM_PARAMETERS = parameters.name != 'null' ? "-P $parameters" : "-p ${sample.getAt(0)}.hmm_params"

    """
    hmm_mrs \\
     -x ${sample.getAt(1)} \\
     -y ${sample.getAt(0)} \\
     -c ${params.MIN_COVERAGE} \\
     -o ${sample.getAt(0)}.MRs.bed \\
     -d ${params.DESERT_SIZE} \\
     -i 30 \\
     -m ${params.MERGE_DIST} \\
     -t ${params.TRIM_METHRATE/100} \\
     ${HUMAN} \\
     ${MIN_C} \\
     ${matrixWG} \\
     ${HMM_PARAMETERS}

     echo -e "${sample.getAt(0)}\t" \\
        \$(cat ${sample.getAt(0)}.MRs.bed | wc -l)"\t" \\
        \$(awk -v OFS=\"\t\" '{sum+=\$3-\$2+1}END{print sum, sprintf("%.0f", sum/NR)}' ${sample.getAt(0)}.MRs.bed) \\
        > ${sample.getAt(0)}.MR_stats.tsv
    """
}

process MethylScore_igv {
    tag "$bed"
    publishDir "${params.PROJECT_FOLDER}/igv", mode: 'copy'

    input:
    file(bed) from MRs_igv.collect()
    file(matrixWG) from matrixWG_igv

    output:
    file('methinfo.igv') into igv

    when:
    params.IGV
 
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
    file(samplesheet) from indexedSamples_splitting.unique().collectFile(){ record -> [ 'samples.tsv', record[0] + '\t' + (record[1]+3) + '\t' + record[0] + '.MRs.bed' + '\n' ] }
    file(MRfile) from MRs_splitting.collect()

    output:
    file('MRbatch*') into MRchunks mode flatten
    file('samples.tsv') into (samples_callDMRs, samples_mergeDMRs)

    script:
    """
    split_MRfile $samplesheet MRbatch ${params.MR_BATCH_SIZE}
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
    file('DMRs.bed') into DMRs
    file('all_context_DMRs.bed') into all_context_DMRs

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
     ${params.HDMR_FOLD_CHANGE}

    sort -k1,1V -k2 -o DMRs.bed DMRs.bed
    sort -k1,1V -k2 -o all_context_DMRs.bed all_context_DMRs.bed
    """
}

workflow.onComplete {
    println()
    if ( workflow.success ) {
     println "[$workflow.complete] >> MethylScore finished SUCCESSFULLY after $workflow.duration and found ${DMRs.getVal().countLines()} DMRs"
    } else {
     println "[$workflow.complete] >> MethylScore finished with ERRORS after $workflow.duration"
    }
}
