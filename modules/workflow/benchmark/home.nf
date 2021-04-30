#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CONVERT  } from '../../process/benchmark/bedGraph2format'
include { RUN_HOME } from '../../process/benchmark/HOME_pairwise'

workflow HOME {
    take:
    samples

    main:

    def contexts = Channel.fromList(params.DMR_CONTEXTS.tokenize(','))

    CONVERT(
        samples,
        "HOME"
    )
    
    RUN_HOME(
       CONVERT.out.converted.groupTuple(by:[0]),
       CONVERT.out.converted.groupTuple(by:[0,1]).collectFile(cache:true, newLine:true){ chrom, group, sample, context, path -> ['samples.tsv', [ group, path.join('\t') ].join('\t') ] },
       contexts
    )
}
