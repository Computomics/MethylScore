workflow GET_SAMPLES {
    main:

    Channel
        .fromPath(params.GENOME, checkIfExists: true)
        .splitFasta( record: [id: true, text: true] )
        .collectFile(storeDir: "${workDir}/fasta") { fasta -> ["${fasta.id}.fa", fasta.text] }
        .set { fasta }

    Channel
        .fromPath(params.SAMPLE_SHEET, checkIfExists: true)
        .splitText()
        .map{ line ->
                def list = line.split()
                def path = file(list[-1])
                assert path.exists(), "${path} doesn't exist!"
                [ 'sampleID':list[0], 'filePath':path ]
            }
        .set { inputMap }

    inputMap
        .map { record -> [ record.sampleID, record.filePath ] }
        .set {samples}

    emit:
    samples
    fasta
}