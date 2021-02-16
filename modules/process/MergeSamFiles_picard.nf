process MERGE_BAM {
    tag "${sampleID}:${bam}"
    publishDir "${params.PROJECT_FOLDER}/01mappings/${sampleID}", mode: 'copy'

    input:
    tuple val(sampleID), path(bam)

    output:
    tuple val(sampleID), path('*.merged.bam'), emit: bam
    
    when:
    bam.size() > 1

    script:
    def mem = task.memory ? task.memory.toMega() : (4.GB).toMega()
    """
    picard \\
        -Xmx${mem}m \\
        -Xms${mem / 4}m \\
        MergeSamFiles \\
            ASSUME_SORTED=true \\
            I=${bam.join(' I=')} \\
            O=${sampleID}.merged.bam
    """
}