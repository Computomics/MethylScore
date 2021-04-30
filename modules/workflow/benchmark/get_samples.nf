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
        //samplesheet needs 4 columns: group, sampleID, context, filepath
        .map{ line -> [ line[0], line[1], line[2], file(line[3], checkIfExists:true) ] }
        //.map{ line -> [ "${line[0]}":line[1], "ctx":line[2], "bedgraph":file(line[3], checkIfExists:true) ] }
        .set { samples }

    emit:
    fasta
    samples
}