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

// validate parameters
ParameterChecks.checkParams(params)

log.info"""
====================================================================================================================================
███╗   ███╗███████╗████████╗██╗  ██╗██╗   ██╗██╗     ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
████╗ ████║██╔════╝╚══██╔══╝██║  ██║╚██╗ ██╔╝██║     ██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
██╔████╔██║█████╗     ██║   ███████║ ╚████╔╝ ██║     ███████╗██║     ██║   ██║██████╔╝█████╗  
██║╚██╔╝██║██╔══╝     ██║   ██╔══██║  ╚██╔╝  ██║     ╚════██║██║     ██║   ██║██╔══██╗██╔══╝  
██║ ╚═╝ ██║███████╗   ██║   ██║  ██║   ██║   ███████╗███████║╚██████╗╚██████╔╝██║  ██║███████╗
╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝   ╚═╝   ╚══════╝╚══════╝ ╚═════╝ ╚═════╝ ╚═╝  ╚═╝╚══════╝v${workflow.manifest.version}
====================================================================================================================================
Reference genome          : ${params.GENOME}
Current home              : $HOME
Current user              : $USER
Current path              : $PWD
Script dir                : $baseDir
Working dir               : $workDir
Output dir                : ${params.PROJECT_FOLDER}
------------------------------------------------------------------------------------------------------------------------------------
SAMPLE_SHEET              : ${params.SAMPLE_SHEET}
------------------------------------------------------------------------------------------------------------------------------------
DO_DEDUP                  : ${params.DO_DEDUP}
HUMAN                     : ${params.HUMAN}
IGV OUTPUT                : ${params.IGV}
ROI                       : ${params.ROI}
MR_PARAMS                 : ${params.MR_PARAMS}
STATISTICS                : ${params.STATISTICS}
------------------------------------------------------------------------------------------------------------------------------------
DMRS_PER_CONTEXT          : ${params.DMRS_PER_CONTEXT}
DMR_CONTEXTS              : ${params.DMRS_CONTEXTS}
CLUSTER_MIN_METH          : ${params.CLUSTER_MIN_METH}
CLUSTER_MIN_METH_CG       : ${params.CLUSTER_MIN_METH_CG}
CLUSTER_MIN_METH_CHG      : ${params.CLUSTER_MIN_METH_CHG}
CLUSTER_MIN_METH_CHH      : ${params.CLUSTER_MIN_METH_CHH}
CLUSTER_MIN_METH_DIFF     : ${params.CLUSTER_MIN_METH_DIFF}
CLUSTER_MIN_METH_DIFF_CG  : ${params.CLUSTER_MIN_METH_DIFF_CG}
CLUSTER_MIN_METH_DIFF_CHG : ${params.CLUSTER_MIN_METH_DIFF_CHG}
CLUSTER_MIN_METH_DIFF_CHH : ${params.CLUSTER_MIN_METH_DIFF_CHH}
DESERT_SIZE               : ${params.DESERT_SIZE}
DMR_MIN_C                 : ${params.DMR_MIN_C}
DMR_MIN_COV               : ${params.DMR_MIN_COV}
FDR_CUTOFF                : ${params.FDR_CUTOFF}
HDMR_FOLD_CHANGE          : ${params.HDMR_FOLD_CHANGE}
IGNORE_FIRST_BP           : ${params.IGNORE_FIRST_BP}
IGNORE_LAST_BP            : ${params.IGNORE_LAST_BP}
MERGE_DIST                : ${params.MERGE_DIST}
MIN_COVERAGE              : ${params.MIN_COVERAGE}
MIN_QUAL                  : ${params.MIN_QUAL}
MR_BATCH_SIZE             : ${params.MR_BATCH_SIZE}
MR_FREQ_CHANGE            : ${params.MR_FREQ_CHANGE}
MR_FREQ_DISTANCE          : ${params.MR_FREQ_DISTANCE}
MR_MIN_C                  : ${params.MR_MIN_C}
SLIDING_WINDOW_SIZE       : ${params.SLIDING_WINDOW_SIZE}
SLIDING_WINDOW_STEP       : ${params.SLIDING_WINDOW_STEP}
TRIM_METHRATE             : ${params.TRIM_METHRATE}
------------------------------------------------------------------------------------------------------------------------------------
Config Profile : ${workflow.profile}
====================================================================================================================================
""".stripIndent()

roi_file = params.ROI ? Channel.fromPath(params.ROI, checkIfExists: true).collect() : file('null')

hmm_params_file = params.MR_PARAMS ? Channel.fromPath(params.MR_PARAMS, checkIfExists: true).collect() : file('null')

/*
 * Create a channel for the reference genome and split it by chromosome
 */

Channel
 .fromPath(params.GENOME, checkIfExists: true)
 .splitFasta( record: [id: true, text: true] )
 .set { fasta }

(fasta_MethylExtract, fasta_bedGraph) = !params.BEDGRAPH ? [ fasta, Channel.empty() ] : [ Channel.empty(), fasta ]

/*
 * Create a Channel for the tsv file containing samples
 * Store entries in a map first, then subset to arrays
 */

Channel
  .fromPath(params.SAMPLE_SHEET, checkIfExists: true)
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
    read_stats.sh \\
     ${sampleID} \\
     ${bamFile} \\
     ${REGIONS_FILE}

    cov_stats.sh \\
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
    paste <(printf "${sampleID.join('\n')}") \\
          <(printf "${consensus.join('\n')}") \\
          <(printf "${Index.join('\n')}") | sort -k 3 > ${chromosomeID}_samples.tsv

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
    python matrix2igv.py -i ${matrixWG} -m MRs.merged.bed -o methinfo.igv
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
     -B $baseDir/bin/betabin_model \\
     -Y $baseDir/bin/pv2qv.py \\
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
     $baseDir/bin/pv2qv.py \\
     ${params.FDR_CUTOFF} \\
     ${params.CLUSTER_MIN_METH} \\
     ${params.DMR_MIN_C} \\
     ${params.HDMR_FOLD_CHANGE}

    sort -k1,1V -k2 -o DMRs.bed DMRs.bed
    sort -k1,1V -k2 -o all_context_DMRs.bed all_context_DMRs.bed
    """
}

workflow.onComplete {
 
    if ( workflow.success ) {
      log.info "[$workflow.complete] >> MethylScore finished SUCCESSFULLY after $workflow.duration and found ${DMRs.getVal().countLines()} DMRs"
    } else {
      log.info "[$workflow.complete] >> MethylScore finished with ERRORS after $workflow.duration"
    }
}
