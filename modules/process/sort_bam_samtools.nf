process SORT {
    tag "$sampleID"
    publishDir "${params.PROJECT_FOLDER}/01mappings/${sampleID}", mode: 'copy'

    input:
    tuple val(sampleID), path(bam)

    output:
    tuple val(sampleID), path("${sampleID}.sorted.bam"), emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sampleID}.sorted.bam $bam
    """
}  