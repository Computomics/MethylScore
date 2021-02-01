process METHYLDACKEL {
    tag "$sampleID:$fasta.id"
    publishDir "${params.PROJECT_FOLDER}/02consensus/${sampleID}/${fasta.id}", mode: 'copy'

    input:
    tuple val(fasta), val(sampleID), path(bamsplit)

    output:
    tuple val(fasta), val(sampleID), path('*.allC'), emit: consensus

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

    sort -k2,2n <(tail -n+2 -q *bedGraph) | awk '{print "${sampleID}\\t"\$0}' > ${sampleID}.allC
    """
}
