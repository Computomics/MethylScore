process SPLIT_MRS {
    tag "batchsize:${params.MR_BATCH_SIZE}"
    publishDir "${params.PROJECT_FOLDER}/05DMRs/batches", mode: 'copy'

    input:
    path(bed)
    path(samplesheet)

    output:
    path('*MRbatch*'), emit: chunks

    script:
    comp = params.PAIRWISE ? "${samplesheet.name.minus('.tsv')}.MRbatch" : "all.MRbatch"
    """
    perl ${projectDir}/bin/split_MRfile.pl ${samplesheet} "${comp}" ${params.MR_BATCH_SIZE}
    """
}