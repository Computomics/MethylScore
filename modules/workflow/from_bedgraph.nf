include { GET_SAMPLES } from './get_samples'
include { MATRIX      } from './get_genome_matrix'
include { SPLIT       } from '../process/split_bedgraph_by_chromosome'

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
    matrixWG = MATRIX.out.matrixWG
    matrixCHROM = MATRIX.out.matrixCHROM
}