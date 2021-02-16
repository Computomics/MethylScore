include { GET_SAMPLES } from './sub/get_samples'
include { SPLIT       } from '../process/split_bedgraph_by_chromosome'
include { MATRIX      } from './sub/get_genome_matrix'

workflow BEDGRAPH {
    main:

    GET_SAMPLES()

    GET_SAMPLES.out.samples
        .groupTuple(by:0)
        .set { samples }

    SPLIT(
        samples,
        GET_SAMPLES.out.fasta
    )

    MATRIX(
        SPLIT.out.consensus,
        samples
    )

    emit:
    matrixCHROM = MATRIX.out.matrixCHROM
    matrixINDEX = MATRIX.out.matrixINDEX
}