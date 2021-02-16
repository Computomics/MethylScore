
include { SPLIT } from '../process/split_matrix_by_chromosome'

workflow MATRIX {   
    main:

    Channel
        .fromPath(params.MATRIX, checkIfExists:true)
        .set { matrix }

    SPLIT(matrix)

    SPLIT.out.matrixCHROM
     .flatten()
     .map { matrix -> [ matrix.name.minus('.genome_matrix.tsv'), matrix ] }
     .set { matrixCHROM }

    emit:
    matrixCHROM
}