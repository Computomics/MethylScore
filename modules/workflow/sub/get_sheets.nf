workflow GENERATE_SHEETS {
    take:
    samples

    main:

    samples
        .flatten()
        .map { it -> def record = it.name.toString().tokenize('__'); [ record[0], record[1] as int ] }
        .set { indexedSamples }

    if ( !params.PAIRWISE ) {

        indexedSamples
            .collectFile(newLine:true,sort:'index'){ sample, index -> ['samples.tsv', [ sample, (index+3), sample + '.MRs.bed' ].join('\t')] }
            .set { mrsheet }

        indexedSamples
            .collectFile(newLine:true,sort:'index'){ sample, index -> ['samples.tsv', [ sample, (index+3) ].join('\t')] }
            .set { dmrsheet }

    } else {

        indexedSamples
            .branch{ sample, index ->
                compare: sample==params.PAIRWISE
                others: sample!=params.PAIRWISE
            }
            .set { branched }

        pairwise = branched.compare.combine(branched.others)

        pairwise
            .collectFile(newLine:true, sort:'index'){ pair -> ["${pair[0]}<>${pair[2]}.tsv", [ [pair[0], pair[1], pair[0] + '.MRs.bed' ].join('\t'), [pair[2], pair[3], pair[2] + '.MRs.bed'].join('\t')].join('\n')] }
            .set { mrsheet }

        pairwise
            .collectFile(newLine:true,sort:'index'){ pair -> ["${pair[0]}<>${pair[2]}.tsv", [ [pair[0], pair[1]].join('\t'), [pair[2], pair[3]].join('\t') ].join('\n')] }
            .set { dmrsheet }
    }

    emit:
    indexedSamples
    mrsheet
    dmrsheet
}