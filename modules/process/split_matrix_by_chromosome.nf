process SPLIT {
    tag "$matrix"
    publishDir "${params.PROJECT_FOLDER}/03matrix/", mode: 'copy'

    input:
    path(matrix)

    output:
    path("*.genome_matrix.tsv"), emit: matrixCHROM

    script:
    """
    awk 'FNR==1{header=\$0;next}{if (end != \$1) close(end); print header > \$1".genome_matrix.tsv"; print >> \$1".genome_matrix.tsv"; end = \$1}' ${matrix}
    """
}
