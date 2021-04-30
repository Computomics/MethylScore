#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { GET_SAMPLES } from './modules/workflow/benchmark/get_samples'
include { METHYLSCORE } from './modules/workflow/benchmark/methylscore'
include { METHYLPY    } from './modules/workflow/benchmark/methylpy'
include { METHYLKIT   } from './modules/workflow/benchmark/methylkit'
include { METILENE    } from './modules/workflow/benchmark/metilene'
include { DSS         } from './modules/workflow/benchmark/dss'
include { DMRSEQ      } from './modules/workflow/benchmark/dmrseq'
include { HOME        } from './modules/workflow/benchmark/home'

workflow {
    
    GET_SAMPLES()

    METHYLSCORE(
        GET_SAMPLES.out.samples,
        GET_SAMPLES.out.fasta
    )
    
    GET_SAMPLES.out.samples.combine(GET_SAMPLES.out.fasta) | (METILENE & DMRSEQ & DSS & METHYLKIT)
    GET_SAMPLES.out.samples.groupTuple(by:[0,1]).combine(GET_SAMPLES.out.fasta) | (METHYLPY & HOME)
}