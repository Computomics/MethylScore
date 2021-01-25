include { INDEX           } from '../process/build_sample_index'
include { GENERATE_SHEETS } from './sub/get_sheets'

workflow MATRIX {   
    main:

    Channel
        .fromPath(params.MATRIX)
        .set { matrixWG }

    INDEX(matrixWG)

    GENERATE_SHEETS(INDEX.out.indexedSamples)

    emit:
    matrixWG
    indexedSamples = GENERATE_SHEETS.out.indexedSamples
    mrsheet = GENERATE_SHEETS.out.mrsheet
    dmrsheet = GENERATE_SHEETS.out.dmrsheet
    index = INDEX.out.index
}
