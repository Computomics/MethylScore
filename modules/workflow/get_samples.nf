workflow GET_SAMPLES {
    main:

    Channel
        .fromPath(params.GENOME, checkIfExists: true)
        .splitFasta( record: [id: true, text: true] )
        .collectFile(cache:true, storeDir: "${workDir}/fasta"){ fasta -> ["${fasta.id}.fa", fasta.text] }
        .map{ fasta -> ["id":fasta.baseName, "seq":fasta] }
        .set { fasta }

    Channel
        .fromPath(params.SAMPLE_SHEET, checkIfExists: true)
        .splitCsv(header:false, strip:true, sep:'\t')
        .map{ line -> [ line[0], file(line[-1], checkIfExists:true) ] }
        .set { samples }

    emit:
    fasta
    samples
}