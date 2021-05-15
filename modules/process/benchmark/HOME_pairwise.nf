process RUN_HOME {
    tag "${fasta.id}"
    publishDir "${params.PROJECT_FOLDER}/HOME/dmrs", mode: 'copy'
    container "quay.io/beckerlab/home-container:latest"

    // symlink data from container into workDir
    beforeScript 'for dir in {saved_model,scripts,training_data};do ln -s /opt/HOME/${dir} ${dir};done'

    input:
    tuple val(fasta), val(group), val(sampleID), val(ctx), path(inputfile)
    path(samplesheet)
    each context

    output:
    path("results_${context}/**/**/*.txt"), emit: home_dmrs

    script:
    """
    HOME-pairwise \\
        -npp ${task.cpus} \\
        -d ${params.CLUSTER_MIN_METH_DIFF/100} \\
        -mc ${params.DMR_MIN_C} \\
        -sin \\
        -t ${context} \\
        -i samples.tsv \\
        -o ./results_${context}
    """
}
