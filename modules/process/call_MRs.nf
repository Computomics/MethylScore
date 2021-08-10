process CALL_MRS {
    tag "${sampleID}:${matrix}"
    publishDir "${params.PROJECT_FOLDER}/04MRs", mode: 'copy', enabled: !params.MR_PARAMS, saveAs: { filename -> filename.endsWith(".hmm_params") ? "hmm_parameters/${filename}" : null }

    input:
    tuple val(sampleID), val(sampleIDX), val(chromID), path(matrix)
    path(parameters)

    output:
    tuple val(chromID), val(sampleID), path("*${sampleID}.MRs.bed"), emit: bed
    path("${sampleID}.hmm_params"), optional: true, emit: hmm_params

    script:
    def HUMAN = params.HUMAN ? "-human" : ""
    def MIN_C = params.MR_MIN_C > 0 ? "-n ${params.MR_MIN_C}" : "-n -1"
    def HMM_PARAMETERS = params.MR_PARAMS ? "-P $parameters" : "-p ${sampleID}.hmm_params"
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
        ${matrix} \\
        ${HMM_PARAMETERS}
    """
}