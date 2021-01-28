process SPLIT {
    tag "$matrix"
    publishDir "${params.PROJECT_FOLDER}/03matrix/", mode: 'copy'

    input:
    path(matrix)

    output:
    path("*.genome_matrix.tsv"), emit: matrixCHROM

    script:
    """
    awk '\$0!~/^#/{if (end != \$1) close(end); print >> \$1".genome_matrix.tsv"; end = \$1}' ${matrix}
    """
}
