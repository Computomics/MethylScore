process SPLIT_MRS {
    tag "batchsize:${params.MR_BATCH_SIZE}"
    publishDir "${params.PROJECT_FOLDER}/05DMRs/batches", mode: 'copy'

    input:
    path(bed)
    path(samplesheet)

    output:
    tuple val(comp), path('*.MRbatch.*'), emit: chunks

    script:
    comp = params.PAIRWISE ? "${samplesheet.name.minus('.tsv')}" : "all"
    """
    for i in *bed; do sort -k1,1 -k2,2n \$i -o \$i; done
    
    perl ${projectDir}/bin/split_MRfile.pl ${samplesheet} "${comp}.MRbatch" ${params.MR_BATCH_SIZE}
    """
}