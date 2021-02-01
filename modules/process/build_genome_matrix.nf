process BUILD {
    tag "$chromosomeID"
    publishDir "${params.PROJECT_FOLDER}/03matrix", mode: 'copy'
    cache 'deep'

    input:
    tuple val(chromosomeID), val(sampleID), path(consensus), path(fasta)
    path(samplesheet)

    output:
    path("*.genome_matrix.tsv"), emit: matrix

    script:
    def input_format = params.METHYLPY ? "methylpy" : "bedgraph"
    """
    sort --parallel=${task.cpus} -m -k3,3g -T . ${consensus} > ${chromosomeID}.allC

    python ${projectDir}/bin/generate_genome_matrix.py \\
        --samples ${samplesheet} \\
        --allC ${chromosomeID}.allC \\
        --format ${input_format} \\
        --fasta ${chromosomeID}.fa
    """
}