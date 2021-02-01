workflow GENERATE_SHEETS {
    take:
    samples

    main:

    samples
        .flatten()
        .map { sample -> def record = sample.name.toString().tokenize('__'); [ record[0], record[1] as int ] }
        .set { indexedSamples }

    if ( !params.PAIRWISE ) {

        indexedSamples
            .collectFile(cache:true, newLine:true, sort:'index'){ sample, index -> ['all.tsv', [ sample, (index+3), sample + '.MRs.bed' ].join('\t')] }
            .set { mrsheet }

        indexedSamples
            .collectFile(cache:true, newLine:true, sort:'index'){ sample, index -> ['all.tsv', [ sample, (index+3) ].join('\t')] }
            .set { dmrsheet }

    } else {

        indexedSamples
            .branch{ sample, index ->
                compare: sample==params.PAIRWISE.toString()
                others:  sample!=params.PAIRWISE.toString()
            }
            .set { branched }

        pairwise = branched.compare.combine(branched.others)

        pairwise
            .collectFile(cache:true, newLine:true, sort:'index'){ pair -> ["${pair[0]}.x.${pair[2]}.tsv", [ [pair[0], (pair[1]+3), pair[0] + '.MRs.bed' ].join('\t'), [pair[2], (pair[3]+3), pair[2] + '.MRs.bed'].join('\t')].join('\n')] }
            .set { mrsheet }

        pairwise
            .collectFile(cache:true, newLine:true, sort:'index'){ pair -> ["${pair[0]}.x.${pair[2]}.tsv", [ [pair[0], (pair[1]+3)].join('\t'), [pair[2], (pair[3]+3)].join('\t') ].join('\n')] }
            .set { dmrsheet }

    }

    emit:
    indexedSamples
    mrsheet
    dmrsheet
}