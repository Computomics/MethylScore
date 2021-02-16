process GENERATE {
    input:
    path(matrix)

    output:
    path('*__*'), emit: indexedSamples

    script:
    """
    INDEX=1
    for i in \$(head -n1 ${matrix} | cut -f 5-); do
        printf "\${i}\\t\${INDEX}\\n" > \${i}__\${INDEX};
        let INDEX++;
    done
    """
}