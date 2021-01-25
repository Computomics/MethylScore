process MERGE_DMRS {
    tag "$comp:$context"
    publishDir "${params.PROJECT_FOLDER}/05DMRs", mode: 'copy'

    input:
    tuple val(comp), path(segments)
    path(samplesheet)
    each context

    output:
    path('*.bed'), emit: dmrs

    script:
    """
    for segment in *.${context}.out/segments.dif; do
        if [ -f \$segment ]; then
            cat \$segment >> ${comp}.${context}.dif;
        else
            touch ${comp}.${context}.dif;
        fi
    done

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

    perl ${projectDir}/bin/merge_DMRs-contexts.pl \\
     ${samplesheet} \\
     ${comp}.${context}.dif \\
     . \\
     parameters.config \\
     ${context}

    sort -k1,1V -k2,2n -o DMRs.${context}.bed DMRs.${context}.bed
    """
}