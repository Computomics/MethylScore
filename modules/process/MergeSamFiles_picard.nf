process MERGE_BAM {
    tag "${sampleID}:${bam}"
    label "resource_medium"
    publishDir "${params.PROJECT_FOLDER}/01mappings/merged/${sampleID}", mode: 'copy', enabled: !params.REMOVE_INTMED_FILES

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