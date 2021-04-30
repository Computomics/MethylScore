#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CONVERT       } from '../../process/benchmark/bedGraph2format'
include { RUN_METHYLKIT } from '../../process/benchmark/methylKit_DMRs'

workflow METHYLKIT {
    take:
    samples

    main:
    CONVERT(
        samples,
        "methylKit"
    )

    RUN_METHYLKIT(
        CONVERT.out.converted.groupTuple(by:[0,3])
    )
}
