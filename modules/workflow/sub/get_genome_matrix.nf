include { BUILD           } from '../../process/build_genome_matrix'
include { INDEX           } from '../../process/build_sample_index'
include { GENERATE_SHEETS } from './get_sheets'

workflow MATRIX {
    take:
    consensus
    fasta
    samples

    main:

    samples
        .collectFile(newLine:true, sort:'index'){ sample, bam -> ['samples.tsv', [ sample, sample + '.allC' ].join('\t')] }
        .set { matrixsheet }

    BUILD(
        consensus.groupTuple(by:0).join(fasta, by:0),
        matrixsheet.collect()
    )

    BUILD.out.matrix
        .collectFile(name: 'genome_matrix.tsv', keepHeader: true, sort: { it.baseName }, storeDir: "${params.PROJECT_FOLDER}/03matrix")
        .set{ matrixWG }

    INDEX(matrixWG)

    GENERATE_SHEETS(INDEX.out.indexedSamples)

    emit:
    matrixWG
    matrixCHROM = BUILD.out.matrix
    indexedSamples = GENERATE_SHEETS.out.indexedSamples
    mrsheet = GENERATE_SHEETS.out.mrsheet
    dmrsheet = GENERATE_SHEETS.out.dmrsheet
    index = INDEX.out.index
}