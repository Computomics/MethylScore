process RUN_DMRFIND {
    tag "${fasta.id}"
    publishDir "${params.PROJECT_FOLDER}/methylpy/dmrs", mode: 'copy'
    container "quay.io/biocontainers/methylpy:1.4.3--py37h41a55b7_0"

    input:
    tuple val(fasta), val(group), val(sampleID), val(ctx), path(allc)
    each context

    output:
    path("*.DMR.bed"), emit: methylpy_dmrs

    script:
    """
    methylpy DMRfind \\
        --allc-files ${allc.join(' ')} \\
        --samples ${sampleID.join(' ')} \\
        --sample-category ${group.join(' ')} \\
        --mc-type ${context} \\
        --chroms ${fasta.id} \\
        --sig-cutoff ${params.FDR_CUTOFF} \\
        --num-procs ${task.cpus} \\
        --min-num-dms ${params.DMR_MIN_C} \\
        --output-prefix methylpy_${context}
    """
}