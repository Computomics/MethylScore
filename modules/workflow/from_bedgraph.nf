include { GET_SAMPLES    } from './get_samples'
include { MATRIX         } from './get_genome_matrix'
include { SPLIT_BEDGRAPH } from '../process/split_bedgraph_by_chromosome'

workflow BEDGRAPH {
    main:

    GET_SAMPLES()

    GET_SAMPLES.out.samples
        .groupTuple(by:0)
        .set { samples }

    SPLIT_BEDGRAPH(
        samples,
        GET_SAMPLES.out.fasta
    )

    MATRIX(
        SPLIT_BEDGRAPH.out.consensus,
        samples
    )

    emit:
    matrixWG = MATRIX.out.matrixWG
    matrixCHROM = MATRIX.out.matrixCHROM
}