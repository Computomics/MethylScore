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

// WIP : nextflow method to split reference by chromosome
//chrfa = Channel.from(reference)
//	       .splitFasta( record: [id: true, seqString: true] )
//	       .collect()
//	       .subscribe{ record -> println record.id }

process MethylScore_deduplicate {
    tag "$sampleID"
    publishDir "${params.PROJECT_FOLDER}/01mappings", mode: 'copy'

    input:
    file reference
    set val(sampleID), val(seqType), file(bamFile) from samples

    output:
    set val(sampleID), val(seqType) into flagW mode flatten
    set val(sampleID), file('*/split/*/*bam'), file('*/split/*/*fa') into fileSplit mode flatten
    file '*/*bam' into dedupBam
    file '*/*bai' into dedupBamIndex
    file '*/*tsv' into readStats
    file '*/*txt' into dedupMetrics

    script:
    """
    mkdir extbin
    ln -s \$(which samtools) extbin/st
    ln -s \$(which bamtools) extbin/bt
    ln -s \$(which bedtools) extbin/bdt
    ln -s /usr/local/share/picard-2.18.11-0/picard.jar extbin/picard.jar
    
    ${baseDir}/${params.SCRIPT_PATH}/merge_and_dedup.sh \\
		1 \\
		. \\
		${sampleID} \\
		${bamFile.join(',')} \\
		${reference} \\
		./extbin \\
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

//get chromosomeID
fileSplit
 .map { chrSplit ->
	chrSplit << chrSplit[1].baseName.split("_")[-1] }
 .set { chrSplit }

process MethylScore_callConsensus {
    tag "$bam"
    publishDir "${params.PROJECT_FOLDER}/02consensus", mode: 'copy'

    input:
    set val(sampleID), file(bam), file(ref), val(chromosome), val(seqType) from chrSplit.combine(flagW.unique(), by: 0)

    output:
    set val(sampleID), file("*/*/${sampleID}.${chromosome}.allC.output"), val(chromosome) into consensus

    script:
// TODO: currently, the reference genome is split in each iteration which could probably be avoided by channeling pre-split chromosomes together with bams using a .join() call
// TODO: chromosome names should be stored in a indexed hashmap, to streamline the sorting steps below
    """
    mkdir extbin
    ln -s \$(which MethylExtract.pl) extbin/MethylExtract.pl

    mkdir --parents ${sampleID}/${chromosome}

    ${baseDir}/${params.SCRIPT_PATH}/consensus.sh \\
		2 \\
		./extbin \\
		${sampleID} \\
		${seqType} \\
		${params.MIN_QUAL} \\
		${params.IGNORE_LAST_BP} \\
		${params.IGNORE_FIRST_BP} \\
		. \\
		${sampleID}/${chromosome} \\
		${params.REMOVE_INTMED_FILES}

    mv ${sampleID}/${chromosome}/allC.output ${sampleID}/${chromosome}/${sampleID}.${chromosome}.allC.output
    """
}

process MethylScore_chromosomalmatrix {
    tag "$chromosome"
    publishDir "${params.PROJECT_FOLDER}/03matrix", mode: 'copy'

    input:
    set val(sampleID), file(pile), val(chromosome) from consensus.groupTuple(by: 2, sort: 'true')

    output:
    set val(chromosome), file("genome_matrix.${chromosome}.tsv") into splitMatrix 
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

splitMatrix
 .transpose()
 .toSortedList({a, b -> a[0] <=> b[0]})
 .flatten()
 .buffer(size:1, skip:1)
 .set { matrix }

idx
 .flatten()
 .unique()
 .toSortedList().withIndex().set{ indexedSamples }

process MethylScore_genomematrix {
    tag "$matrixlist"
    publishDir "${params.PROJECT_FOLDER}/03matrix", mode: 'copy'

    input:
    file(matrixlist) from matrix.collect()

    output:
    file('genome_matrix.tsv') into matrixWG

    script:
    """
    cat ${matrixlist} | sed '1!{/^#/d;}' > genome_matrix.tsv
    """
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

    ${baseDir}/${params.SCRIPT_PATH}/call_DMRs.sh 5 '${baseDir}/${params.BIN_PATH}/dmrs -s ${samples} -r ${chunk} -m ${matrix} -p ${params.MR_FREQ_CHANGE} -i ${params.CLUSTER_MIN_METH_DIFF} -j ${params.CLUSTER_MIN_METH} -v ${params.DMR_MIN_COV} -n ${params.DMR_MIN_C} -w ${params.SLIDING_WINDOW_SIZE} -x ${params.SLIDING_WINDOW_STEP} -z 1 -B ${baseDir}/${params.BIN_PATH}/betabin_model -Y ${baseDir}/${params.SCRIPT_PATH}/pv2qv.py --no-post-process -o ${chunk}.out' ${chunk}.out
    """
}

process MethylScore_mergeDMRs {
    tag "$segments"
    publishDir "${params.PROJECT_FOLDER}/05DMRs", mode: 'copy'

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
