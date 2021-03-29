process MERGE_DMRS {
    tag "${context}:${samplesheet.baseName}"
    publishDir "${params.PROJECT_FOLDER}/05DMRs/${samplesheet.baseName}", mode: 'copy'

    input:
    tuple val(comp), val(context), path(segments), path(samplesheet)

    output:
    path("*.bed"), emit: dmrs

    script:
    """
    cat <<EOF >> parameters.config
    PYTHON_PATH:  "python"
    SCRIPT_PATH: "${projectDir}/bin"
    FDR_CUTOFF: ${params.FDR_CUTOFF}
    CLUSTER_MIN_METH_CG:  ${params.CLUSTER_MIN_METH_CG}
    CLUSTER_MIN_METH_CHG: ${params.CLUSTER_MIN_METH_CHG}
    CLUSTER_MIN_METH_CHH: ${params.CLUSTER_MIN_METH_CHH}
    DMR_MIN_C:  ${params.DMR_MIN_C}
    HDMR_FOLD_CHANGE: ${params.HDMR_FOLD_CHANGE}
    EOF

    merge_DMRs-contexts.pl \\
        ${samplesheet} \\
        ${segments} \\
        . \\
        parameters.config \\
        ${context}

    sort -k1,1V -k2,2n -o DMRs.${context}.bed DMRs.${context}.bed
    """
}