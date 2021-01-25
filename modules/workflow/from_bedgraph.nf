include { GET_SAMPLES } from './sub/get_samples'
include { SPLIT       } from '../process/split_bedgraph_by_chromosome'
include { INDEX       }  from '../process/build_sample_index'
include { MATRIX      } from './sub/get_genome_matrix'

workflow BEDGRAPH {
    main:

    GET_SAMPLES()

    SPLIT(
        GET_SAMPLES.out.samples,
        GET_SAMPLES.out.fasta
    )

    MATRIX(
        SPLIT.out.consensus,
        SPLIT.out.fasta,
        GET_SAMPLES.out.samples
    )

    emit:
    matrixWG = MATRIX.out.matrixWG
    indexedSamples = MATRIX.out.indexedSamples
    mrsheet = MATRIX.out.mrsheet
    dmrsheet = MATRIX.out.dmrsheet
    index = MATRIX.out.index
}