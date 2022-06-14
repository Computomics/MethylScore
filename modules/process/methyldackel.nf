process METHYLDACKEL {
    tag "${sampleID}:${fasta.id}"
    label "resource_low"
    publishDir "${params.PROJECT_FOLDER}/02consensus", mode: 'copy', saveAs: { filename -> filename.endsWith(".svg") ? "mbias/${filename}" : !params.REMOVE_INTMED_FILES ? "${sampleID}/${fasta.id}/${filename}" : null }

    input:
    tuple val(fasta), val(sampleID), path(bamsplit)

    output:
    tuple val(fasta), val(sampleID), path('*.allC'), emit: consensus
    path('*.svg'), optional: true, emit: mbias

    script:
    """
    MethylDackel extract \\
        --CHH \\
        --CHG \\
        --nOT ${params.IGNORE_OT} --nOB ${params.IGNORE_OB} \\
        -p ${params.MIN_QUAL} \\
        --minOppositeDepth=1 \\
        --maxVariantFrac=0.01 \\
        --keepDupes \\
        ${fasta.seq} \\
        ${bamsplit}

    MethylDackel mbias \\
        --CHH \\
        --CHG \\
        ${fasta.seq} \\
        ${bamsplit} \\
        ${sampleID}.${fasta.id}

    sort -k2,2n <(tail -n+2 -q *bedGraph) | awk '{print "${sampleID}\\t"\$0}' > ${sampleID}.allC
    """
}
