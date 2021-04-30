process CONVERT {
    tag "${context}"
    publishDir "${params.PROJECT_FOLDER}/${format}", mode: 'copy'
    
    input:
    tuple val(group), val(sampleID), val(context), path(bedGraph), val(fasta)
    val(format)

    output:
    tuple val(fasta), val(group), val(sampleID), val(context), path("*"), emit: converted

    script:
    def ctx  = (context.toList().size > 1) ? "combined" : context
    def inputfile = (bedGraph.toList().size > 1) ? "<(sort -m -k2,2n ${bedGraph})" : "${bedGraph}"

    if (format == 'metilene')
    """
    cut -f-4 ${bedGraph} > ${bedGraph.baseName}.${format}.${ctx}.input
    """
    else
    """
    python ${projectDir}/bin/convert_bedgraph.py \\
        --format ${format} \\
        --input ${inputfile} \\
        --output ${sampleID}.${format}.${ctx}.input \\
        --sample ${sampleID} \\
        --fasta ${fasta.seq}
    """
}