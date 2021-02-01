process SPLIT {
    tag "$sampleID:$fasta.id"
    publishDir "${params.PROJECT_FOLDER}/01mappings/${sampleID}/split/${fasta.id}/", mode: 'copy'

    input:
    tuple val(sampleID), path(bedgraph)
    each fasta

    output:
    tuple val(fasta), val(sampleID), path('*.allC'), optional: true, emit: consensus

    script:
    def compressed = bedgraph.first().toString().endsWith('gz') ? 'zcat' : 'cat'
    """
    ${compressed} ${bedgraph} | awk '\$1 == "${fasta.id}"' | sort -k2,2n | awk '{ if(NR > 0) {print "${sampleID}\\t"\$0 > "${sampleID}.allC"}}'
    """
}
