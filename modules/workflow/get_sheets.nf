include { GENERATE } from '../process/get_sample_index'

workflow SAMPLESHEET {
    take:
    matrixWG

    main:
    GENERATE(matrixWG)

    GENERATE.out.indexedSamples
        .flatten()
        .map { sample -> def record = sample.name.toString().tokenize('__'); [ record[0], record[1] as int ] }
        .set { indexedSamples }

    if ( !params.PAIRWISE ) {

        indexedSamples
            .collectFile(cache:true, newLine:true, sort:'index'){ sample, index -> ["all.tsv", [ sample, (index+3), sample + '.MRs.bed' ].join('\t')] }
            .set { sheet }

    } else {

        indexedSamples
            .branch{ sample, index ->
                A: sample==params.PAIRWISE.toString()
                B:  sample!=params.PAIRWISE.toString()
            }
            .set { branched }

        pairwise = branched.A.combine(branched.B)

        pairwise
            .collectFile(cache:true, newLine:true, sort:'index'){ pair -> ["${pair[0]}.x.${pair[2]}.tsv", [ [pair[0], (pair[1]+3), pair[0] + '.MRs.bed' ].join('\t'), [pair[2], (pair[3]+3), pair[2] + '.MRs.bed'].join('\t')].join('\n')] }
            .set { sheet }

    }

    emit:
    indexedSamples
    sheet
}