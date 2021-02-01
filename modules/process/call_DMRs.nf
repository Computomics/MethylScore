process CALL_DMRS {
    tag "$context:$chunk"
    publishDir "${params.PROJECT_FOLDER}/05DMRs/${comp}/batches", mode: 'copy'
    cache 'deep'

    input:
    tuple path(matrixWG), path(tbi)
    path(samplesheet)
    tuple val(comp), path(chunk)
    each context

    output:
    tuple val(comp), val(context), path("${chunk}.${context}.out/*dif"), optional: true, emit: segments

    script:
    def CLUSTER_MIN_METH = !params.DMRS_PER_CONTEXT ? params.CLUSTER_MIN_METH : params."CLUSTER_MIN_METH_${context}"
    def CLUSTER_MIN_METH_DIFF = !params.DMRS_PER_CONTEXT ? params.CLUSTER_MIN_METH_DIFF : params."CLUSTER_MIN_METH_DIFF_${context}"

    """
    perl ${projectDir}/bin/SEGMENTS-contexts.pl \\
     -c ${context} \\
     -s ${comp}.tsv \\
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