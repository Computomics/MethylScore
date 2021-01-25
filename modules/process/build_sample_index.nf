process INDEX {
    tag "${matrixWG}"
    publishDir "${params.PROJECT_FOLDER}/03matrix", mode: 'copy'

    input:
    path(matrixWG)

    output:
    tuple path('*.gz'), path('*.tbi'), emit: index
    path('*__*'), emit: indexedSamples

    script:
    """
    INDEX=1
        for i in \$(head -n1 ${matrixWG} | cut -f 5-); do
        printf "\${i}\\t\${INDEX}\\n" > \${i}__\${INDEX};
        let INDEX++;
    done

    bgzip -f --threads ${task.cpus} ${matrixWG}
    tabix -s 1 -b 2 -e 2 ${matrixWG}.gz
    """
}