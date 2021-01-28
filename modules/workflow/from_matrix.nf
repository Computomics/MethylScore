include { INDEX           } from '../process/build_sample_index'
include { SPLIT           } from '../process/split_matrix_by_chromosome'
include { GENERATE_SHEETS } from './sub/get_sheets'

workflow MATRIX {   
    main:

    Channel
        .fromPath(params.MATRIX)
        .set { matrixWG }

    matrixWG | (INDEX & SPLIT)

    GENERATE_SHEETS(INDEX.out.indexedSamples)

    emit:
    matrixWG
    matrixCHROM = SPLIT.out.matrixCHROM.flatten()
    indexedSamples = GENERATE_SHEETS.out.indexedSamples
    mrsheet = GENERATE_SHEETS.out.mrsheet
    dmrsheet = GENERATE_SHEETS.out.dmrsheet
    index = INDEX.out.index
}
