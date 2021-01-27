process SORT {
    tag "$sampleID"
    publishDir "${params.PROJECT_FOLDER}/01mappings/${sampleID}", mode: 'copy'

    input:
    tuple val(sampleID), path(bam)

    output:
    tuple val(sampleID), path("*.sorted.bam"), emit: bam

    script:
    """
    samtools view @ ${task.cpus} -bh -F 0x200 -F 0x4 | samtools sort -@ ${task.cpus} -o ${bam.baseName}.sorted.bam $bam
    """
} 