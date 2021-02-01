process CALL_MRS {
    tag "${sampleID}:${matrix}"
    publishDir "${params.PROJECT_FOLDER}/04MRs/${sampleID}", mode: 'copy'
    cache 'deep'

    input:
    tuple val(sampleID), val(sampleIDX)
    each path(matrix)
    path(parameters)

    output:
    tuple val(sampleID), path("*${sampleID}.MRs.bed"), emit: bed
    path("${sampleID}.hmm_params"), optional: true, emit: hmm_params
    path("*${sampleID}.MR_stats.tsv"), emit: stats

    script:
    def HUMAN = params.HUMAN ? "-human" : ""
    def MIN_C = params.MR_MIN_C > 0 ? "-n ${params.MR_MIN_C}" : "-n -1"
    def HMM_PARAMETERS = params.MR_PARAMS ? "-P $parameters" : "-p ${sampleID}.hmm_params"
    def CHROM = params.MR_PARAMS ? "${matrix.name.toString().tokenize('.')[0]}." : ""
    """
    hmm_mrs \\
     -x ${sampleIDX} \\
     -y ${sampleID} \\
     -c ${params.MIN_COVERAGE} \\
     -o ${CHROM}${sampleID}.MRs.bed \\
     -d ${params.DESERT_SIZE} \\
     -i 30 \\
     -m ${params.MERGE_DIST} \\
     -t ${params.TRIM_METHRATE/100} \\
     ${HUMAN} \\
     ${MIN_C} \\
     ${matrix} \\
     ${HMM_PARAMETERS}

    MR_stats.sh ${CHROM}${sampleID} ${CHROM}${sampleID}.MRs.bed
    """
}