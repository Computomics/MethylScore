
include { SPLIT } from '../process/split_matrix_by_chromosome'

workflow MATRIX {   
    main:

    Channel
        .fromPath(params.MATRIX, checkIfExists:true)
        .map { matrix -> [ 'all', matrix ] }
        .set { matrixWG }

    matrixCHROM = params.MR_PARAMS ? matrixWG | SPLIT | flatten | map { matrix -> [ matrix.name.minus('.genome_matrix.tsv'), matrix ] } : Channel.empty()

    emit:
    matrixWG
    matrixCHROM
}