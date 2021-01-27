process BUILD {
    tag "$chromosomeID"
    publishDir "${params.PROJECT_FOLDER}/03matrix", mode: 'copy'

    input:
    tuple val(chromosomeID), val(sampleID), path(consensus), path(fasta)
    path(samplesheet)

    output:
    path("*.genome_matrix.tsv"), emit: matrix

    script:
    def input_format = params.METHYLPY ? "methylpy" : "bismark"
    """
    sort --parallel=${task.cpus} -m -k3,3g -T . ${consensus} > ${chromosomeID}.allC

    perl ${projectDir}/bin/generate_genome_matrix.pl \\
        -s ${samplesheet} \\
        -f ${chromosomeID}.allC \\
        -i ${input_format} \\
        -r ${chromosomeID}.fa \\
        -o ${chromosomeID}.genome_matrix.tsv
    """
}