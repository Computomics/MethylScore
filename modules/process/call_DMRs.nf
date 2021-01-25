process CALL_DMRS {
    tag "$context:$chunk"
    publishDir "${params.PROJECT_FOLDER}/05DMRs/${comp}", mode: 'copy'

    input:
    tuple path(matrixWG), path(tbi)
    path(samplesheet)
    path(chunk)
    each context

    output:
    tuple val(comp), path("${chunk}.${context}.out"), optional: true, emit: segments

    script:
    comp = params.PAIRWISE ? chunk.toString().tokenize('.MRbatch')[0] : 'all'

    def CLUSTER_MIN_METH = !params.DMRS_PER_CONTEXT ? params.CLUSTER_MIN_METH : params."CLUSTER_MIN_METH_${context}"
    def CLUSTER_MIN_METH_DIFF = !params.DMRS_PER_CONTEXT ? params.CLUSTER_MIN_METH_DIFF : params."CLUSTER_MIN_METH_DIFF_${context}"

    """
    perl ${projectDir}/bin/SEGMENTS-contexts.pl \\
     -c ${context} \\
     -s ${samplesheet} \\
     -r ${chunk} \\
     -m ${matrixWG} \\
     -p ${params.MR_FREQ_CHANGE} \\
     -i ${CLUSTER_MIN_METH_DIFF} \\
     -j ${CLUSTER_MIN_METH} \\
     -v ${params.DMR_MIN_COV} \\
     -n ${params.DMR_MIN_C} \\
     -w ${params.SLIDING_WINDOW_SIZE} \\
     -x ${params.SLIDING_WINDOW_STEP} \\
     -z 1 \\
     -B $projectDir/bin/betabin_model \\
     -Y $projectDir/bin/pv2qv.py \\
     --no-post-process \\
     -o ${chunk}.${context}.out
    """
}