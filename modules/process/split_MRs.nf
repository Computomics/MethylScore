process SPLIT_MRS {
    tag "${chromID}:batchsize:${params.MR_BATCH_SIZE}"
    label "resource_medium"
    publishDir "${params.PROJECT_FOLDER}/05DMRs/${comp}/batches", mode: 'copy', enabled: !params.REMOVE_INTMED_FILES

    input:
    tuple val(chromID), val(sampleID), path(bed)
    each path(samplesheet)

    output:
    tuple val(chromID), val(comp), path('*.MRbatch.*'), path(samplesheet), emit: chunks

    script:
    comp = samplesheet.name.minus('.tsv')
    """
    split_MRfile.pl \\
        ${samplesheet} \\
        "${comp}.MRbatch" \\
        ${params.MR_BATCH_SIZE}
    """
}