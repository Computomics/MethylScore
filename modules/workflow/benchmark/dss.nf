#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CONVERT } from '../../process/benchmark/bedGraph2format'
include { RUN_DSS } from '../../process/benchmark/dss_DMRs'

workflow DSS {
    take:
    samples

    main:
    CONVERT(
        samples,
        "DSS"
    )

    RUN_DSS(
        CONVERT.out.converted.groupTuple(by:[0,3])
    )
}
