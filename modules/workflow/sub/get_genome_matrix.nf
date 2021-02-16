include { BUILD } from '../../process/build_genome_matrix'

workflow MATRIX {
    take:
    consensus
    samples

    main:

    samples
        .collectFile(cache:true, newLine:true, sort:'index'){ sample, bam -> ['samples.tsv', [ sample, sample + '.allC' ].join('\t')] }
        .set { matrixsheet }

    BUILD(
        consensus.groupTuple(by:0, sort:true),
        matrixsheet.collect()
    )

    emit:
    matrixCHROM = BUILD.out.matrix
}