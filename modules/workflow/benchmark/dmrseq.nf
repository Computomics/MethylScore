#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { RUN_DMRSEQ } from '../../process/benchmark/dmrseq_DMRs'

workflow DMRSEQ {
    take:
    samples

    main:

    RUN_DMRSEQ(
        samples.groupTuple(by:2)
    )
}
