process SPLIT_MRS {
    tag "${chromID}:batchsize:${params.MR_BATCH_SIZE}"
    publishDir "${params.PROJECT_FOLDER}/05DMRs/batches", mode: 'copy'

    input:
    tuple val(chromID), val(sampleID), path(bed)
    each path(samplesheet)

    output:
    tuple val(chromID), val(comp), path('*.MRbatch.*'), path(samplesheet), emit: chunks

    script:
    comp = samplesheet.name.minus('.tsv')
    """
    perl ${projectDir}/bin/split_MRfile.pl \\
        ${samplesheet} \\
        "${comp}.MRbatch" \\
        ${params.MR_BATCH_SIZE}
    """
}