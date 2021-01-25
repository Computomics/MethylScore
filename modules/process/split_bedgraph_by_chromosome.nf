process SPLIT {
    tag "$sampleID:$chromosomeID"
    publishDir "${params.PROJECT_FOLDER}/01mappings/${sampleID}/split/${chromosomeID}/", mode: 'copy'

    input:
    tuple val(sampleID), path(bedgraph)
    each path(fasta)

    output:
    tuple val(chromosomeID), val(sampleID), path('*.allC'), emit: consensus
    tuple val(chromosomeID), path(fasta), emit: fasta

    script:
    chromosomeID = fasta.baseName
    def compressed = bedgraph.first().toString().endsWith('gz') ? 'zcat' : 'cat'
    """
    ${compressed} ${bedgraph} | awk '\$1 == "${chromosomeID}"' | sort -k2,2n | awk '{print "${sampleID}\\t"\$0}' > ${sampleID}.allC
    """
}
