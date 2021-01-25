process CALL_MRS {
    tag "${sampleIDX}:${sampleID}"
    publishDir "${params.PROJECT_FOLDER}/04MRs/${sampleID}", mode: 'copy'

    input:
    tuple val(sampleID), val(sampleIDX)
    path(matrixWG)
    path(parameters)

    output:
    tuple val(sampleID), path("${sampleID}.MRs.bed"), emit: bed
    path("${sampleID}.hmm_params"), optional: true, emit: hmm_params
    path("${sampleID}.MR_stats.tsv"), emit: stats

    script:
    def HUMAN = params.HUMAN ? "-human" : ""
    def MIN_C = params.MR_MIN_C > 0 ? "-n ${params.MR_MIN_C}" : "-n -1"
    def HMM_PARAMETERS = parameters.name != 'null' ? "-P $parameters" : "-p ${sampleID}.hmm_params"

    """
    hmm_mrs \\
     -x ${sampleIDX} \\
     -y ${sampleID} \\
     -c ${params.MIN_COVERAGE} \\
     -o ${sampleID}.MRs.bed \\
     -d ${params.DESERT_SIZE} \\
     -i 30 \\
     -m ${params.MERGE_DIST} \\
     -t ${params.TRIM_METHRATE/100} \\
     ${HUMAN} \\
     ${MIN_C} \\
     ${matrixWG} \\
     ${HMM_PARAMETERS}

    MR_stats.sh ${sampleID} ${sampleID}.MRs.bed
    """
}