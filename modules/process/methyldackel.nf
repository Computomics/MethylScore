process METHYLDACKEL {
    tag "$sampleID:$chromosomeID"
    publishDir "${params.PROJECT_FOLDER}/02consensus/${sampleID}/${chromosomeID}", mode: 'copy'

    input:
    tuple val(chromosomeID), val(sampleID), path(bamsplit), path(fasta)

    output:
    tuple val(chromosomeID), val(sampleID), path('*.allC'), emit: consensus

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
    ${chromosomeID}.fa \\
    ${bamsplit}

    sort -k2,2n <(tail -n+2 -q *bedGraph) | awk '{print "${sampleID}\\t"\$0}' > ${sampleID}.allC
    """
}
