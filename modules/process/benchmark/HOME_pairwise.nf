process RUN_HOME {
    tag "${fasta.id}"
    publishDir "${params.PROJECT_FOLDER}/HOME/dmrs", mode: 'copy'
    container '/scratch-cbe/users/patrick.huether/simulated/home_b6796e84.sif'

    // symlink data from container into workDir
    beforeScript 'for dir in {saved_model,scripts,training_data};do ln -s /opt/HOME/${dir} ${dir};done'
    //beforeScript 'ln -s /opt/HOME/saved_model saved_model && ln -s /opt/HOME/scripts scripts && ln -s /opt/HOME/training_data training_data'

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
