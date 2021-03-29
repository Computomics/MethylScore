include { BUILD } from '../process/build_genome_matrix'

workflow MATRIX {
    take:
    consensus
    samples

    main:

    samples
        .collectFile(cache:true, newLine:true, sort:'index'){ sample, consensus -> ['samples.tsv', [ sample, sample + '.allC' ].join('\t')] }
        .set { matrixsheet }

    BUILD(
        consensus.groupTuple(by:0, sort:true),
        matrixsheet.collect()
    )

    BUILD.out.matrix
        .collectFile(cache:true, keepHeader:true, sort:{ it.baseName }, storeDir:"${params.PROJECT_FOLDER}/03matrix"){ chromID, matrix -> [ 'genome_matrix.tsv', matrix ]}
        .map { matrix -> [ 'all', matrix ] }
        .set{ matrixWG }

    emit:
    matrixWG
    matrixCHROM = BUILD.out.matrix
}