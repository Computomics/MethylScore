process SPLIT {
    tag "$sampleID:$chromosomeID"
    publishDir "${params.PROJECT_FOLDER}/01mappings/${sampleID}/split/${chromosomeID}/", mode: 'copy'

    input:
    tuple val(sampleID), path(bam)
    each path(fasta)

    output:
    tuple val(chromosomeID), val(sampleID), path("${chromosomeID}.bam"), emit: bam
    tuple val(chromosomeID), path(fasta), emit: fasta

    script:
    chromosomeID = fasta.baseName
    """
    samtools index ${bam}
    cat <(samtools view -H ${bam} | grep -E "@HD|SN:${chromosomeID}\$(printf '\\t')") \\
        <(samtools view ${bam} ${chromosomeID}) \\
        | samtools view -bo ${chromosomeID}.bam -
    """
}
