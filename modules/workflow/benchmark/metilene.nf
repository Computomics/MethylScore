#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CONVERT      } from '../../process/benchmark/bedGraph2format'
include { RUN_METILENE } from '../../process/benchmark/metilene_DMRs'

workflow METILENE {
    take:
    samples

    main:
    CONVERT(
        samples,
        "metilene"
    )

    RUN_METILENE(
        CONVERT.out.converted.groupTuple(by:[0,3])
    )
}
