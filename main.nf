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

===================================================================================================================================
Reference genome      : ${params.GENOME}
Current home          : $HOME
Current user          : $USER
Current path          : $PWD
Script dir            : $baseDir
Working dir           : $workDir
Output dir            : ${params.PROJECT_FOLDER}
-----------------------------------------------------------------------------------------------------------------------------------
SAMPLE_SHEET          : ${params.SAMPLE_SHEET}
-----------------------------------------------------------------------------------------------------------------------------------
DO_DEDUP              : ${params.DO_DEDUP}
HUMAN                 : ${params.HUMAN}
IGV OUTPUT            : ${params.IGV}
ROI                   : ${params.ROI}
MR_PARAMS             : ${params.MR_PARAMS}
STATISTICS            : ${params.STATISTICS}
-----------------------------------------------------------------------------------------------------------------------------------
CLUSTER_MIN_METH      : ${params.CLUSTER_MIN_METH}
CLUSTER_MIN_METH_DIFF : ${params.CLUSTER_MIN_METH_DIFF}
DESERT_SIZE           : ${params.DESERT_SIZE}
DMR_MIN_C             : ${params.DMR_MIN_C}
DMR_MIN_COV           : ${params.DMR_MIN_COV}
FDR_CUTOFF            : ${params.FDR_CUTOFF}
HDMR_FOLD_CHANGE      : ${params.HDMR_FOLD_CHANGE}
IGNORE_OT             : ${params.AUTOTRIM ? 'auto' : params.IGNORE_OT}
IGNORE_OB             : ${params.AUTOTRIM ? 'auto' : params.IGNORE_OB}
MERGE_DIST            : ${params.MERGE_DIST}
MIN_COVERAGE          : ${params.MIN_COVERAGE}
MIN_QUAL              : ${params.MIN_QUAL}
MR_BATCH_SIZE         : ${params.MR_BATCH_SIZE}
MR_FREQ_CHANGE        : ${params.MR_FREQ_CHANGE}
MR_FREQ_DISTANCE      : ${params.MR_FREQ_DISTANCE}
MR_MIN_C              : ${params.MR_MIN_C}
SLIDING_WINDOW_SIZE   : ${params.SLIDING_WINDOW_SIZE}
SLIDING_WINDOW_STEP   : ${params.SLIDING_WINDOW_STEP}
TRIM_METHRATE         : ${params.TRIM_METHRATE}
-----------------------------------------------------------------------------------------------------------------------------------
Config Profile : ${workflow.profile}
===================================================================================================================================
""".stripIndent()

if( params.AUTOTRIM ){
log.warn "MethylScore is running in AUTOTRIM mode. Please review mbias plots in ${params.PROJECT_FOLDER}/01mappings/mbias and adjust IGNORE_OT and IGNORE_OB settings if necessary"
}

roi_file = params.ROI ? Channel.fromPath(params.ROI, checkIfExists: true).collect() : file('null')
hmm_params_file = params.MR_PARAMS ? Channel.fromPath(params.MR_PARAMS, checkIfExists: true).collect() : file('null')

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
        [ 'sampleID':list[0], 'filePath':file(list[-1]) ]
      }
  .set { inputMap }

def sampleIndex = 1
inputMap
  .map { record -> [ record.sampleID, record.filePath ] }
  .groupTuple()
  .tap { samples }
  .map { record -> tuple( record[0], sampleIndex++ ) }
  .set{ indexedSamples }

(samples_bam, samples_bedGraph) = !params.BEDGRAPH ? [ samples, Channel.empty() ] : [ Channel.empty(), samples ]

/*
 * Start pipeline
 */

process MethylScore_mergeReplicates {
    tag "$sampleID: ${bamFile.collect().size()} replicate(s)"
    publishDir "${params.PROJECT_FOLDER}/01mappings/${sampleID}", mode: 'copy'

    input:
    set val(sampleID), file(bamFile) from samples_bam

    output:
    set val(sampleID), file('*.bam') into mergedSamples

    when:
    !params.BEDGRAPH

    script:
    if( bamFile.toList().size() != 1 )
       """
       picard -Xmx${task.memory.toMega()}m -Xms${task.memory.toMega() / 4}m -Djava.io.tmpdir='.' -XX:ParallelGCThreads=1 \\
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
    set val(sampleID), file(bamFile) from mergedSamples

    output:
    set val(sampleID), file('*dedup.bam') into (dedup, read_stats)
    file ('dedup.metrics.txt')

    when:
    !params.BEDGRAPH

    script:
    """
    picard -Xmx${task.memory.toMega()}m -Xms${task.memory.toMega() / 4}m -Djava.io.tmpdir='.' -XX:ParallelGCThreads=1 \\
      MarkDuplicates \\
        I=${bamFile} \\
        O=${sampleID}.dedup.bam \\
        METRICS_FILE=dedup.metrics.txt \\
        REMOVE_DUPLICATES=true \\
        MAX_FILE_HANDLES=1 \\
        TMP_DIR=. \\
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
    set val(sampleID), file(bamFile) from read_stats
    file(ROIs) from roi_file

    output:
    file ('*') into stats

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
    publishDir "${params.PROJECT_FOLDER}/01mappings", mode: 'copy',
       saveAs: {filename -> filename.endsWith(".svg") ? "mbias/$filename" : "${sampleID}/split/${chromosomeID}/${filename}"}

    input:
    set val(sampleID), file(bamFile) from dedup.mix(samples_bedGraph)
    each file(fasta) from fasta_split

    output:
    set val(sampleID), file("${sampleID}.${chromosomeID}.{bam,allC}"), val(chromosomeID) into chrSplit
    set val(sampleID), stdout, val(chromosomeID) optional true into mbias
    file('*.svg') optional true into mbias_plots

    script:
    chromosomeID = fasta.baseName

    if (!params.BEDGRAPH)
      """
      samtools index ${bamFile}
      cat <(samtools view -H ${bamFile} | grep -E '@HD|${chromosomeID}') \\
          <(samtools view ${bamFile} ${chromosomeID}) | \\
          samtools view -bo ${sampleID}.${chromosomeID}.bam -

      MethylDackel mbias \\
          --CHH \\
          --CHG \\
          ${chromosomeID}.fa \\
          ${sampleID}.${chromosomeID}.bam \\
          ${sampleID}.${chromosomeID} 2> >(tail -n1 | cut -d: -f2) > /dev/null
      """
    else
      """
      awk '\$1 == "${chromosomeID}"' ${bamFile} | sort -k2,2n > ${sampleID}.${chromosomeID}.allC
      """
}

(bamSplit, bedSplit) = !params.BEDGRAPH ? [ chrSplit, Channel.empty() ] : [ Channel.empty(), chrSplit ]

process MethylScore_callConsensus {
    tag "$sampleID:$chromosomeID"
    publishDir "${params.PROJECT_FOLDER}/02consensus/${sampleID}/${chromosomeID}", mode: 'copy'

    input:
    set val(sampleID), val(chromosomeID), file(splitBam), val(mbias) from bamSplit.combine(mbias, by:[0,2])
    file(fasta) from fasta_consensus.collect()

    output:
    set val(sampleID), file('*.allC'), val(chromosomeID) into allC

    when:
    !params.BEDGRAPH

    script:
    def trim = params.AUTOTRIM ? "${mbias.trim()}" : "--nOT ${params.IGNORE_OT} --nOB ${params.IGNORE_OB}" 

    """
    MethylDackel extract \\
     --CHH \\
     --CHG \\
     ${trim} \\
     -p ${params.MIN_QUAL} \\
     --minOppositeDepth=1 \\
     --maxVariantFrac=0.01 \\
     --keepDupes \\
     ${chromosomeID}.fa \\
     ${splitBam}

    tail -n+2 -q *bedGraph | sort -k2,2n > ${sampleID}.${chromosomeID}.allC
    """
}

consensus = !params.BEDGRAPH ? allC : bedSplit

indexedSamples
 .combine(consensus, by: 0)
 .groupTuple(by: 3)
 .set {pile}

process MethylScore_chromosomalmatrix {
    tag "${chromosomeID}"
    publishDir "${params.PROJECT_FOLDER}/03matrix", mode: 'copy'

    input:
    set val(sampleID), val(Index), file(consensus), val(chromosomeID) from pile
    file(fasta) from fasta_matrix.collect()

    output:
    file("${chromosomeID}.genome_matrix.tsv") into splitMatrix
    set val(sampleID), val(Index) into (indexedSamples_splitting, indexedSamples_MRs) mode flatten

    script:

    """
    paste <(printf "${sampleID.join('\n')}") \\
          <(printf "${consensus.join('\n')}") \\
          <(printf "${Index.join('\n')}") | sort -k 3,3n > ${chromosomeID}_samples.tsv

    generate_genome_matrix \\
     -s ${chromosomeID}_samples.tsv \\
     -i bismark -r ${chromosomeID}.fa \\
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
    file("${sample.getAt(0)}.MR_stats.tsv") into MR_stats

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
    file(bed) from MRs_igv.collect()
    file(matrixWG) from matrixWG_igv

    output:
    file('methinfo.igv') into igv

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
    file('*/*.dif') into segmentFiles

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

    sort -k1,1V -k2,2n -o DMRs.bed DMRs.bed
    sort -k1,1V -k2,2n -o all_context_DMRs.bed all_context_DMRs.bed
    """
}
