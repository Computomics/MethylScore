process SPLIT_MATRIX {
    tag "$matrix"
    publishDir "${params.PROJECT_FOLDER}/03matrix/", mode: 'copy'

    input:
    tuple val(chromID), path(matrix)

    output:
    path("*.genome_matrix.tsv"), emit: matrixCHROM

    shell:
    '''
    awk -F '\t' 'FNR==1{header=$0;next}{if (!seen[$1]++) print header > $1".genome_matrix.tsv"; print >> $1".genome_matrix.tsv"}' !{matrix}
    '''
}
