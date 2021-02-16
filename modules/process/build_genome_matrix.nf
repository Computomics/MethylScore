process BUILD {
    tag "${fasta.id}"
    publishDir "${params.PROJECT_FOLDER}/03matrix/", mode: 'copy', enabled: !params.REMOVE_INTMED_FILES

    input:
    tuple val(fasta), val(sampleID), path(consensus)
    path(samplesheet)

    output:
    tuple val(fasta.id), path("*.genome_matrix.tsv"), emit: matrix

    script:
    def input_format = params.METHYLPY ? "methylpy" : "bedgraph"
    """
    sort --parallel=${task.cpus} -m -k3,3g -T . ${consensus} > ${fasta.id}.allC

    python ${projectDir}/bin/generate_genome_matrix.py \\
        --samples ${samplesheet} \\
        --allC ${fasta.id}.allC \\
        --format ${input_format} \\
        --fasta ${fasta.seq}
    """
}