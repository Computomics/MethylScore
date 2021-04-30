#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CONVERT     } from '../../process/benchmark/bedGraph2format'
include { RUN_DMRFIND } from '../../process/benchmark/methylpy_DMRfind'

workflow METHYLPY {
    take:
    samples

    main:

    def contexts = Channel.fromList(params.DMR_CONTEXTS.tokenize(','))

    CONVERT(
        samples,
        "methylpy"
    )

    RUN_DMRFIND(
        CONVERT.out.converted.groupTuple(by:0),
        contexts
    )
}
