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
params.IGV = false
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
                  .groupTuple()

samples.into { samples_mergeAndDedup; samples_consensus }

process MethylScore_deduplicate {
    tag "$sampleID"
    publishDir "${params.PROJECT_FOLDER}/01mappings", mode: 'copy'
    module 'BamTools/2.4.0-foss-2016a:SAMtools/1.3.1-foss-2016a:Java/1.8.0_112'
// TODO: switch to container solution to resolve dependencies

    input:
    file reference
    set val(sampleID), val(seqType), file(bamFile) from samples_mergeAndDedup

    output:
    set val(sampleID), val(seqType), file('*/split/*/*bam') into bamSplit mode flatten
    set val(sampleID), file('*/split/*/*fa') into refSplit mode flatten
    file '*/*bam' into dedupBam
    file '*/*bai' into dedupBamIndex
    file '*/*tsv' into readStats
    file '*/*txt' into dedupMetrics

    script:
    """
    ${baseDir}/${params.SCRIPT_PATH}/merge_and_dedup.sh \\
		1 \\
		. \\
		${sampleID} \\
		${bamFile.join(',')} \\
		${reference} \\
		${baseDir}/${params.EXTBIN_PATH} \\
		${params.FORCE_RERUN} \\
		${params.REMOVE_INTMED_FILES} \\
		1024 \\
		bamtools,split \\
		samtools \\
		${params.STATISTICS} \\
		${baseDir}/${params.SCRIPT_PATH} \\
		${params.DO_DEDUP} \\
		"${params.ROI}"
    """
}

process MethylScore_callConsensus {
    tag "$bam"
    publishDir "${params.PROJECT_FOLDER}/02consensus", mode: 'copy'
    module 'BamTools/2.4.0-foss-2016a:SAMtools/1.3.1-foss-2016a:Perl/5.22.1-foss-2016a'

    input:
    set val(sampleID), val(seqType), file(bam) from bamSplit
    set val(sampleID), file(ref) from refSplit

    output:
    set val(sampleID), file("*/*/${sampleID}.${chromosome}.allC.output"), val(chromosome) into consensus

    script:
    // extract the chromosome ID, as we need it for grouping the tuple lateron
    chromosome = bam.baseName.split("_")[-1]
// TODO: currently, the reference genome is split in each iteration which could probably be avoided by channeling pre-split chromosomes together with bams using a .join() call
// TODO: chromosome names should be stored in a indexed hashmap, to streamline the sorting steps below
    """
    mkdir --parents ${sampleID}/${chromosome}

    ${baseDir}/${params.SCRIPT_PATH}/consensus.sh \\
		2 \\
		${baseDir}/${params.EXTBIN_PATH} \\
		${sampleID} \\
		${seqType} \\
		${params.MIN_QUAL} \\
		${params.IGNORE_LAST_BP} \\
		${params.IGNORE_FIRST_BP} \\
		. \\
		${sampleID}/${chromosome} \\
		${params.REMOVE_INTMED_FILES} \\
		samtools

    mv ${sampleID}/${chromosome}/allC.output ${sampleID}/${chromosome}/${sampleID}.${chromosome}.allC.output
    """
}

process MethylScore_chromosomalmatrix {
    tag "$chromosome"
    publishDir "${params.PROJECT_FOLDER}/03matrix", mode: 'copy'

    input:
    set val(sampleID), file(pile), val(chromosome) from consensus.groupTuple(by: 2, sort: true)

    output:
    file('*tsv') into matrix
    val(sampleID) into idx

    script:
    """
    for i in ${pile}; do
    	echo -e \$(basename \$i | cut -f1 -d.)'\t'\$i >> samples.txt
    done
    sort samples.txt > sorted_samples.txt

    ${baseDir}/${params.SCRIPT_PATH}/generate_matrix.sh 3 ${baseDir}/${params.BIN_PATH} sorted_samples.txt . ${chromosome}
    """
// TODO cutting filenames like this is not robust against samplenames containing dots
// TODO sorting like this is probably unnecessary and correct order of samples can be handled more elegantly in the .collectFile() call below
}

idx
 .flatten()
 .unique()
 .toSortedList().withIndex().set{ indexedSamples }

process MethylScore_genomematrix {
    tag "$matrix"
    publishDir "${params.PROJECT_FOLDER}/03matrix", mode: 'copy'

    input:
    file(matrix) from matrix.collectFile(name:'matrixWG.tsv')

    output:
    file('*tsv') into matrixWG
    file('*header') into matrixheader

    script:
    """
    head -n1 ${matrix} > genome_matrix.header

    sort -k1,1V -k2,2n ${matrix} | sed '/^#/d' > matrix
    cat genome_matrix.header matrix > genome_matrix.tsv
    """
// TODO sorting should be handled in groovy, not bash
}

matrixWG.into{matrixWG_MRs; matrixWG_igv; matrixWG_DMRs}

process MethylScore_callMRs {
    tag "$sample"
    publishDir "${params.PROJECT_FOLDER}/04MRs", mode: 'copy'

    input:
    file(matrixWG) from matrixWG_MRs
    each sample from indexedSamples

    output:
    file("${sample[0]}/${sample[1]+4}.${sample[0]}.MRs.bed") into MRs
    file('*/*.params') into hmmparams
    file('*/*.tsv') into MRstats

    script:
    """
    ${baseDir}/${params.SCRIPT_PATH}/call_MRs.sh \\
		4 \\
		${baseDir}/${params.BIN_PATH} \\
		${sample[0]} \\
		${sample[1]+1} \\
		${params.MIN_COVERAGE} \\
		${params.DESERT_SIZE} \\
		${params.MERGE_DIST} \\
		${params.TRIM_METHRATE/100} \\
		${matrixWG} \\
		./${sample[0]} \\
		${params.HUMAN} \\
		${params.MR_MIN_C}
    mv ${sample[0]}/MRs.bed ${sample[0]}/${sample[1]+4}.${sample[0]}.MRs.bed
    """
}

MRs.into{MRs_igv; MRs_split}

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
  ${baseDir}/${params.SCRIPT_PATH}/igv.sh "" ${baseDir}/${params.SCRIPT_PATH} ${matrixWG} . ${bed}
  """
}

process MethylScore_splitMRs {
    tag "batchsize:${params.MR_BATCH_SIZE}"
    publishDir "${params.PROJECT_FOLDER}/05DMRs/batches", mode: 'copy'

    input:
    file(bed) from MRs_split.collect()

    output:
    file('*') into chunks mode flatten
    file('../samplesheet.tsv') into chunkfile

    script:
    """
    for i in ${bed}; do
      echo -e \$(basename \$i | cut -f2 -d.)'\t'\$(basename \$i | cut -f1 -d.)'\t'\$i >> ../samples.txt
    done
    sort ../samples.txt > ../samplesheet.tsv

    split_MRfile ../samplesheet.tsv MRbatch ${params.MR_BATCH_SIZE}
    """
}

process MethylScore_callDMRs {
    tag "$chunk"
    publishDir "${params.PROJECT_FOLDER}/05DMRs", mode: 'copy'
    module 'Python/2.7.11-foss-2016a'

    input:
    file(matrix) from matrixWG_DMRs
    file(samples) from chunkfile
    each file(chunk) from chunks

    output:
    file('*/*.bed') optional true into beds
    file('*/*.dif') optional true into difs
    file('samplesheet.tsv') into mergefile

    script:
    """
    mkdir ${chunk}.out

    ${baseDir}/${params.SCRIPT_PATH}/call_DMRs.sh 5 '${baseDir}/${params.BIN_PATH}/dmrs -s ${samples} -r ${chunk} -m ${matrix} -p ${params.MR_FREQ_CHANGE} -i ${params.CLUSTER_MIN_METH_DIFF} -j ${params.CLUSTER_MIN_METH} -v ${params.DMR_MIN_COV} -n ${params.DMR_MIN_C} -w ${params.SLIDING_WINDOW_SIZE} -x ${params.SLIDING_WINDOW_STEP} -z 1 -B ${baseDir}/${params.BIN_PATH}/betabin_model -T ${baseDir}/${params.EXTBIN_PATH}/tbx -E ${baseDir}/${params.EXTBIN_PATH}/bzp -K " " -Y ${baseDir}/${params.SCRIPT_PATH}/pv2qv.py --no-post-process -o ${chunk}.out' ${chunk}.out
    """
}

process MethylScore_mergeDMRs {
    tag "$segments"
    publishDir "${params.PROJECT_FOLDER}/05DMRs", mode: 'copy'
    module 'Python/2.7.11-foss-2016a'

    input:
    file(segments) from difs.collectFile(name:'segments.dif')
    file(samples) from mergefile

    output:
    file('*') into dmrs

    script:
    """
    ${baseDir}/${params.BIN_PATH}/merge_DMRs \\
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
