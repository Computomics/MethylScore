workflow GET_SAMPLES {
    main:

    Channel
        .fromPath(params.GENOME, checkIfExists: true)
        .splitFasta( record: [id: true, text: true] )
        .collectFile(storeDir: "${workDir}/fasta") { fasta -> ["${fasta.id}.fa", fasta.text] }
        .set { fasta }

    Channel
        .fromPath(params.SAMPLE_SHEET, checkIfExists: true)
        .splitCsv(header:false, strip:true, sep:'\t')
        .map{ line -> [ line[0], file(line[-1], checkIfExists:true) ] }
        .set { samples }

    emit:
    samples
    fasta
}