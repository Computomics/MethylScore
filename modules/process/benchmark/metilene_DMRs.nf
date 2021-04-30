process RUN_METILENE {
    tag "${fasta.id}"
    publishDir "${params.PROJECT_FOLDER}/metilene/dmrs", mode: 'copy'
    container 'quay.io/biocontainers/metilene:0.2.8--h516909a_0'

    input:
    tuple val(fasta), val(group), val(sampleID), val(context), path(inputfile)

    output:
    path("*.DMRs.bed"), emit: metilene_dmrs

    script:
    // For CHH context, segments tend to be very long,
    // in which case metilene requires excessive amounts of memory.
    // The -G parameter limits the maximum segment size
    def max_length = context == 'CHH' ? '2000' : '-1'
    """
    bedtools unionbedg \\
        -header \\
        -names ${sampleID.join(' ')} \\
        -i ${inputfile.join(' ')} \\
        -filler 'NA' \\
        | cut -f1,3- | sed 's/end/pos/' > metilene.input

    metilene \\
        -t ${task.cpus} \\
        -m ${params.DMR_MIN_C} \\
        -d ${params.CLUSTER_MIN_METH_DIFF} \\
        -a ${group.unique()[0]} \\
        -b ${group.unique()[1]} \\
        -G ${max_length} \\
        metilene.input > metilene.${context}.DMRs.bed
    """
}