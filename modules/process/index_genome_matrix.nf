process BUILD_INDEX {
    tag "$matrix"
    publishDir "${params.PROJECT_FOLDER}/03matrix/", mode: 'copy'

    input:
    tuple val(chromID), path(matrix)

    output:
    tuple val(chromID), path('*.gz'), path('*.tbi'), emit: index

    script:
    """
    bgzip \\
        -f \\
        --threads ${task.cpus} \\
        ${matrix}

    tabix \\
        -s 1 \\
        -b 2 \\
        -e 2 \\
        ${matrix}.gz
    """
}