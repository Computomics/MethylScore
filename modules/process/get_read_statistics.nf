process READ_STATISTICS {
    tag "${sampleID}"
    publishDir "${params.PROJECT_FOLDER}/01mappings/${sampleID}", mode: 'copy', enabled: !params.REMOVE_INTMED_FILES

    input:
    tuple val(sampleID), path(bam)
    path(regions)

    output:
    path('*'), emit: stats

    script:
    def roi = regions.name != 'null' ? "${regions}" : ""

    """
    read_stats.sh ${sampleID} ${bam} ${roi}
    cov_stats.sh ${sampleID} ${bam} ${roi}
    """
}