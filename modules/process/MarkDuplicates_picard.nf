process DEDUPLICATE {
    tag "${sampleID}"
    label "resource_high"
    publishDir "${params.PROJECT_FOLDER}/01mappings/deduplicated/${sampleID}", mode: 'copy', enabled: !params.REMOVE_INTMED_FILES

    input:
    tuple val(sampleID), path(bam)

    output:
    tuple val(sampleID), path('*.dedup.bam') , emit: bam

    script:
    def mem = task.memory ? task.memory.toMega() : (4.GB).toMega()
    """
    picard \\
        -Xmx${mem - 512}m \\
        -Xms${mem / 4}m \\
        -XX:ParallelGCThreads=1 \\
        MarkDuplicates \\
            I=${bam} \\
            O=${sampleID}.dedup.bam \\
            METRICS_FILE=dedup.metrics.txt \\
            REMOVE_DUPLICATES=true \\
            MAX_FILE_HANDLES=1 \\
            VALIDATION_STRINGENCY=LENIENT
    """
}