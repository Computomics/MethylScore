process SORT_BAM {
    tag "${sampleID}"
    label "resource_medium"
    publishDir "${params.PROJECT_FOLDER}/01mappings/sorted/${sampleID}", mode: 'copy', enabled: !params.REMOVE_INTMED_FILES

    input:
    tuple val(sampleID), path(bam)

    output:
    tuple val(sampleID), path("*.sorted.bam"), emit: bam

    script:
    """
    samtools view -@ ${task.cpus} -bh -F 0x200 -F 0x4 $bam | samtools sort -@ ${task.cpus} -o ${bam.baseName}.sorted.bam
    """
}
