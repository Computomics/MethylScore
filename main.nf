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

log.info """
===================================================================================================================================

███╗   ███╗███████╗████████╗██╗  ██╗██╗   ██╗██╗     ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
████╗ ████║██╔════╝╚══██╔══╝██║  ██║╚██╗ ██╔╝██║     ██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
██╔████╔██║█████╗     ██║   ███████║ ╚████╔╝ ██║     ███████╗██║     ██║   ██║██████╔╝█████╗  
██║╚██╔╝██║██╔══╝     ██║   ██╔══██║  ╚██╔╝  ██║     ╚════██║██║     ██║   ██║██╔══██╗██╔══╝  
██║ ╚═╝ ██║███████╗   ██║   ██║  ██║   ██║   ███████╗███████║╚██████╗╚██████╔╝██║  ██║███████╗
╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝   ╚═╝   ╚══════╝╚══════╝ ╚═════╝ ╚═════╝ ╚═╝  ╚═╝╚══════╝v${workflow.manifest.version}

====================================================================================================================================
Current home              : $HOME
Current user              : $USER
Current path              : $PWD
Script dir                : $baseDir
Working dir               : $workDir
------------------------------------------------------------------------------------------------------------------------------------
PROJECT_FOLDER            : ${params.PROJECT_FOLDER}
------------------------------------------------------------------------------------------------------------------------------------
GENOME                    : ${params.GENOME}
SAMPLE_SHEET              : ${params.SAMPLE_SHEET}
------------------------------------------------------------------------------------------------------------------------------------
BEDGRAPH                  : ${params.BEDGRAPH}
DO_DEDUP                  : ${params.BEDGRAPH ? "ignored (BEDGRAPH = ${params.BEDGRAPH})" : params.DO_DEDUP }
STATISTICS                : ${params.BEDGRAPH ? "ignored (BEDGRAPH = ${params.BEDGRAPH})" : params.STATISTICS }
HUMAN                     : ${params.HUMAN}
IGV OUTPUT                : ${params.IGV}
ROI                       : ${params.ROI}
MR_PARAMS                 : ${params.MR_PARAMS}
METRICS                   : ${params.METRICS}
------------------------------------------------------------------------------------------------------------------------------------
DMRS_PER_CONTEXT          : ${params.DMRS_PER_CONTEXT}
DMR_CONTEXTS              : ${params.DMRS_PER_CONTEXT ? params.DMR_CONTEXTS : 'combined'}
CLUSTER_MIN_METH          : ${params.DMRS_PER_CONTEXT ? "ignored (DMRS_PER_CONTEXT = ${params.DMRS_PER_CONTEXT})" : params.CLUSTER_MIN_METH}
CLUSTER_MIN_METH_DIFF     : ${params.DMRS_PER_CONTEXT ? "ignored (DMRS_PER_CONTEXT = ${params.DMRS_PER_CONTEXT})" : params.CLUSTER_MIN_METH_DIFF}
CLUSTER_MIN_METH_CG       : ${!params.DMRS_PER_CONTEXT ? "ignored (DMRS_PER_CONTEXT = ${params.DMRS_PER_CONTEXT})" : params.CLUSTER_MIN_METH_CG}
CLUSTER_MIN_METH_CHG      : ${!params.DMRS_PER_CONTEXT ? "ignored (DMRS_PER_CONTEXT = ${params.DMRS_PER_CONTEXT})" : params.CLUSTER_MIN_METH_CHG}
CLUSTER_MIN_METH_CHH      : ${!params.DMRS_PER_CONTEXT ? "ignored (DMRS_PER_CONTEXT = ${params.DMRS_PER_CONTEXT})" : params.CLUSTER_MIN_METH_CHH}
CLUSTER_MIN_METH_DIFF_CG  : ${!params.DMRS_PER_CONTEXT ? "ignored (DMRS_PER_CONTEXT = ${params.DMRS_PER_CONTEXT})" : params.CLUSTER_MIN_METH_DIFF_CG}
CLUSTER_MIN_METH_DIFF_CHG : ${!params.DMRS_PER_CONTEXT ? "ignored (DMRS_PER_CONTEXT = ${params.DMRS_PER_CONTEXT})" : params.CLUSTER_MIN_METH_DIFF_CHG}
CLUSTER_MIN_METH_DIFF_CHH : ${!params.DMRS_PER_CONTEXT ? "ignored (DMRS_PER_CONTEXT = ${params.DMRS_PER_CONTEXT})" : params.CLUSTER_MIN_METH_DIFF_CHH}
DESERT_SIZE               : ${params.DESERT_SIZE}
DMR_MIN_C                 : ${params.DMR_MIN_C}
DMR_MIN_COV               : ${params.DMR_MIN_COV}
FDR_CUTOFF                : ${params.FDR_CUTOFF}
HDMR_FOLD_CHANGE          : ${params.HDMR_FOLD_CHANGE}
IGNORE_OT                 : ${params.BEDGRAPH ? "ignored (BEDGRAPH = ${params.BEDGRAPH})" : params.IGNORE_OT}
IGNORE_OB                 : ${params.BEDGRAPH ? "ignored (BEDGRAPH = ${params.BEDGRAPH})" : params.IGNORE_OB}
MERGE_DIST                : ${params.MERGE_DIST}
MIN_COVERAGE              : ${params.MIN_COVERAGE}
MIN_QUAL                  : ${params.BEDGRAPH ? "ignored (BEDGRAPH = ${params.BEDGRAPH})" : params.MIN_QUAL }
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

// Create a channel for contexts to be analysed. Set to 'combined' if DMRS_PER_CONTEXT = false, because this is what dmrs-contexts expects
DMRcontexts = params.DMRS_PER_CONTEXT ? Channel.from(params.DMR_CONTEXTS.tokenize(',')) : Channel.from('combined')

/*
 * Create a channel for the reference genome and split it by chromosome
 */

Channel
 .fromPath(params.GENOME, checkIfExists: true)
 .splitFasta( record: [id: true, text: true] )
 .collectFile(storeDir: "${workDir}/fasta") { fasta -> ["${fasta.id}.fa", fasta.text] }
 .into { fasta_split; fasta_consensus; fasta_matrix }

/*
 * Create a Channel for the tsv file containing samples
 * Store entries in a map first, then subset to arrays
 */

Channel
  .fromPath(params.SAMPLE_SHEET, checkIfExists: true)
  .splitText()
  .map{ line ->
        def list = line.split()
        def path = file(list[-1])
        assert path.exists(), "${path} doesn't exist!"
        [ 'sampleID':list[0], 'filePath':path ]
      }
  .set { inputMap }

def sampleIndex = 1
inputMap
  .map { record -> [ record.sampleID, record.filePath ] }
  .groupTuple()
  .tap { samples }
  .map { record -> tuple( record[0], sampleIndex++ ) }
  .set { indexedSamples }

(samples_bam, samples_bedGraph) = !params.BEDGRAPH ? [ samples, Channel.empty() ] : [ Channel.empty(), samples ]

/*
 * Start pipeline
 */

process MethylScore_mergeReplicates {
    tag "$sampleID: ${bamFile.collect().size()} replicate(s)"
    publishDir "${params.PROJECT_FOLDER}/01mappings/${sampleID}", mode: 'copy'

    input:
    tuple val(sampleID), path(bamFile) from samples_bam

    output:
    tuple val(sampleID), path('*.bam') into mergedSamples

    when:
    !params.BEDGRAPH

    script:
    if( bamFile.toList().size() != 1 )
       """
       picard -Xmx${task.memory.toMega()}m -Xms${task.memory.toMega() / 4}m -XX:ParallelGCThreads=1 \\
        MergeSamFiles \\
        I=${bamFile.join(' I=')} \\
        O=${sampleID}.merged.sorted.bam \\
        USE_THREADING=false
       """
    else
       """
       if [[ \$(samtools view -H ${bamFile} | grep unsorted) ]]; then
         samtools sort -o ${sampleID}.sorted.bam ${bamFile};
       else
         mv ${bamFile} ${sampleID}.sorted.bam;
       fi
       """
}

if( params.DO_DEDUP ) {

 process MethylScore_deduplicate {
    tag "$sampleID"
    publishDir "${params.PROJECT_FOLDER}/01mappings/${sampleID}", mode: 'copy'

    input:
    tuple val(sampleID), path(bamFile) from mergedSamples

    output:
    tuple val(sampleID), path('*dedup.bam') into (dedup, read_stats)
    path('dedup.metrics.txt')

    when:
    !params.BEDGRAPH

    script:
    """
    picard -Xmx${task.memory.toMega() - 512}m -Xms${task.memory.toMega() / 4}m -XX:ParallelGCThreads=1 \\
      MarkDuplicates \\
        I=${bamFile} \\
        O=${sampleID}.dedup.bam \\
        METRICS_FILE=dedup.metrics.txt \\
        REMOVE_DUPLICATES=true \\
        MAX_FILE_HANDLES=1 \\
        VALIDATION_STRINGENCY=LENIENT
    """
}

} else {

(dedup, read_stats) = mergedSamples.into(2)

}

process MethylScore_readStatistics {
    tag "$sampleID"
    publishDir "${params.PROJECT_FOLDER}/01mappings/${sampleID}", mode: 'copy'

    input:
    tuple val(sampleID), path(bamFile) from read_stats
    path(ROIs) from roi_file

    output:
    path('*') into stats

    when:
    params.STATISTICS && !params.BEDGRAPH

    script:
    def REGIONS_FILE = ROIs.name != 'null' ? "${ROIs}" : ""

    """
    read_stats.sh ${sampleID} ${bamFile} ${REGIONS_FILE}
    cov_stats.sh ${sampleID} ${bamFile} ${REGIONS_FILE}
    """
}

process MethylScore_splitBams {
    tag "$sampleID:$chromosomeID"
    publishDir "${params.PROJECT_FOLDER}/01mappings/${sampleID}/split/${chromosomeID}/", mode: 'copy'

    input:
    tuple val(sampleID), path(bamFile) from dedup.mix(samples_bedGraph)
    each path(fasta) from fasta_split

    output:
    tuple val(sampleID), path("${sampleID}.{${chromosomeID}.bam,meth}"), val(chromosomeID) into chrSplit

    script:
    chromosomeID = fasta.baseName

    if (!params.BEDGRAPH)
      """
      samtools index ${bamFile}
      cat <(samtools view -H ${bamFile} | grep -E '@HD|${chromosomeID}') \\
          <(samtools view ${bamFile} ${chromosomeID}) | \\
          samtools view -bo ${sampleID}.${chromosomeID}.bam -
      """
    else
      """
      cat ${bamFile} | awk '\$1 == "${chromosomeID}"' > ${sampleID}.${chromosomeID}.meth
      """
}

//(bamSplit, bedSplit) = !params.BEDGRAPH ? [ chrSplit, Channel.empty() ] : [ Channel.empty(), chrSplit ]

process MethylScore_callConsensus {
    tag "$sampleID:$chromosomeID"
    publishDir "${params.PROJECT_FOLDER}/02consensus/${sampleID}/${chromosomeID}", mode: 'copy'

    input:
    tuple val(sampleID), val(chromosomeID), path(splitChrom) from chrSplit
    path(fasta) from fasta_consensus.collect()

    output:
    tuple val(sampleID), path('*.allC'), val(chromosomeID) into consensus

    script:
    if (!params.BEDGRAPH)
      """
      MethylDackel extract \\
      --CHH \\
      --CHG \\
      --nOT ${params.IGNORE_OT} --nOB ${params.IGNORE_OB} \\
      -p ${params.MIN_QUAL} \\
      --minOppositeDepth=1 \\
      --maxVariantFrac=0.01 \\
      --keepDupes \\
      ${chromosomeID}.fa \\
      ${splitChrom}

      LC_ALL=C sort -k2,2n <(tail -n+2 -q *bedGraph) | awk '{print "${sampleID}\\t"\$0}' > ${sampleID}.allC
      """
    else
      """
      LC_ALL=C sort -k2,2n ${splitChrom} | awk '{print "${sampleID}\\t"\$0}' > ${sampleID}.allC
      """   
}

//consensus = !params.BEDGRAPH ? allC : bedSplit

indexedSamples
 .tap { indexedSamples_matrix; indexedSamples_callMRs; indexedSamples_splitMRs; indexedSamples_callDMRs; indexedSamples_mergeDMRs }
 .combine(consensus, by: 0)
 .groupTuple(by: 3)
 .set { pile }

process MethylScore_chromosomalmatrix {
    tag "${chromosomeID}"
    publishDir "${params.PROJECT_FOLDER}/03matrix", mode: 'copy'

    input:
    tuple val(sampleID), val(Index), path(consensus), val(chromosomeID) from pile
    path(fasta) from fasta_matrix.collect()
    path(samples) from indexedSamples_matrix.collectFile(){ record -> [ "samples.tsv", record[0] + '\t' + record[0] + '.allC' + '\n' ] }.collect()

    output:
    path("${chromosomeID}.genome_matrix.tsv") into splitMatrix

    script:
    def input_format = params.METHYLPY ? "methylpy" : "bismark" 
    """
    LC_ALL=C sort -m -T . ${consensus} > merged.tsv
    perl ${baseDir}/bin/generate_genome_matrix.pl \\
     -s ${samples} \\
     -f merged.tsv \\
     -i ${input_format} -r ${chromosomeID}.fa \\
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
    path(matrixWG) from matrixWG_MRs
    each sample from indexedSamples_callMRs
    path(parameters) from hmm_params_file

    output:
    path("${sample.getAt(0)}.MRs.bed") into (MRs_igv, MRs_splitting)
    path("${sample.getAt(0)}.hmm_params") optional true into hmm_params
    path("${sample.getAt(0)}.MR_stats.tsv") into MR_stats

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

    MR_stats.sh ${sample.getAt(0)} ${sample.getAt(0)}.MRs.bed
    """
}

process MethylScore_igv {
    tag "$bed"
    publishDir "${params.PROJECT_FOLDER}/igv", mode: 'copy'

    input:
    path(bed) from MRs_igv.collect()
    path(matrixWG) from matrixWG_igv

    output:
    path('methinfo.igv') into igv

    when:
    params.IGV
 
    script:
    """
    sort -m -k1,1 -k2,2g -k3,3g ${bed} > MRs.merged.bed
    python $baseDir/bin/matrix2igv.py -i ${matrixWG} -m MRs.merged.bed -o methinfo.igv
    """
}

process MethylScore_splitMRs {
    tag "batchsize:${params.MR_BATCH_SIZE}"
    publishDir "${params.PROJECT_FOLDER}/05DMRs/batches", mode: 'copy'

    input:
    path(samples) from indexedSamples_splitMRs.collectFile(){ record -> [ 'samples.tsv', record[0] + '\t' + (record[1]+3) + '\t' + record[0] + '.MRs.bed' + '\n' ] }
    path(MRfile) from MRs_splitting.collect()

    output:
    path('MRbatch*') into MRchunks mode flatten

    script:
    """
    perl ${baseDir}/bin/split_MRfile.pl ${samples} MRbatch ${params.MR_BATCH_SIZE}
    """
}

process MethylScore_callDMRs {
    tag "$context:$chunk"
    publishDir "${params.PROJECT_FOLDER}/05DMRs", mode: 'copy'

    input:
    path(matrixWG) from matrixWG_DMRs
    path(samples) from indexedSamples_callDMRs.collectFile(){ record -> [ 'samples.tsv', record[0] + '\t' + (record[1]+3) + '\n' ] }.collect()
    each path(chunk) from MRchunks
    each context from DMRcontexts

    output:
    path('*/*.bed') optional true into bedFiles
    tuple val(context), path('*/*.dif') optional true into segmentFiles

    script:
    def cluster_min_meth = !params.DMRS_PER_CONTEXT ? params.CLUSTER_MIN_METH : params."CLUSTER_MIN_METH_${context}"
    def cluster_min_meth_diff = !params.DMRS_PER_CONTEXT ? params.CLUSTER_MIN_METH_DIFF : params."CLUSTER_MIN_METH_DIFF_${context}"

    """
    mkdir ${chunk}.${context}.out

    perl ${baseDir}/bin/SEGMENTS-contexts.pl \\
     -c ${context} \\
     -s ${samples} \\
     -r ${chunk} \\
     -m ${matrixWG} \\
     -p ${params.MR_FREQ_CHANGE} \\
     -i ${cluster_min_meth_diff} \\
     -j ${cluster_min_meth} \\
     -v ${params.DMR_MIN_COV} \\
     -n ${params.DMR_MIN_C} \\
     -w ${params.SLIDING_WINDOW_SIZE} \\
     -x ${params.SLIDING_WINDOW_STEP} \\
     -z 1 \\
     -B $baseDir/bin/betabin_model \\
     -Y $baseDir/bin/pv2qv.py \\
     --no-post-process \\
     -o ${chunk}.${context}.out
    """
}

process MethylScore_mergeDMRs {
    tag "$segments.name"
    publishDir "${params.PROJECT_FOLDER}/05DMRs", mode: 'copy'

    input:
    path(segments) from segmentFiles.collectFile(name: { it[0] })
    path(samples) from indexedSamples_mergeDMRs.collectFile(){ record -> [ 'samples.tsv', record[0] + '\t' + (record[1]+3) + '\n' ] }.collect()

    output:
    path('*.bed') into DMRs

    script:
    def context = segments.name
    """
    cat <<EOF >> parameters.config
    PYTHON_PATH:  "python"
    SCRIPT_PATH: "${baseDir}/bin"
    FDR_CUTOFF: ${params.FDR_CUTOFF}
    CLUSTER_MIN_METH_CG:  ${params.CLUSTER_MIN_METH_CG}
    CLUSTER_MIN_METH_CHG: ${params.CLUSTER_MIN_METH_CHG}
    CLUSTER_MIN_METH_CHH: ${params.CLUSTER_MIN_METH_CHH}
    DMR_MIN_C:  ${params.DMR_MIN_C}
    HDMR_FOLD_CHANGE: ${params.HDMR_FOLD_CHANGE}
    EOF

    perl ${baseDir}/bin/merge_DMRs-contexts.pl \\
     ${samples} \\
     ${segments} \\
     . \\
     parameters.config \\
     ${context}

    sort -k1,1V -k2,2n -o DMRs.${context}.bed DMRs.${context}.bed
    """
}

workflow.onComplete {
    if ( workflow.success ) {
      log.info "[$workflow.complete] >> MethylScore finished SUCCESSFULLY after $workflow.duration and found ${DMRs.collectFile(name: 'allDMRs.bed').getVal().countLines()} DMRs"
    } else {
      log.info "[$workflow.complete] >> MethylScore finished with ERRORS after $workflow.duration"
    }
}