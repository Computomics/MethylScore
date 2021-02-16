process SPLIT_BAM {
    tag "${sampleID}:${fasta.id}"
    publishDir "${params.PROJECT_FOLDER}/01mappings/${sampleID}/split/${fasta.id}/", mode: 'copy', enabled: !params.REMOVE_INTMED_FILES

    input:
    tuple val(sampleID), path(bam)
    each fasta

    output:
    tuple val(fasta), val(sampleID), path("${fasta.id}.bam"), emit: bam

    script:
    """
    samtools index ${bam}
    cat <(samtools view -H ${bam} | grep -E "@HD|SN:${fasta.id}\$(printf '\\t')") \\
        <(samtools view ${bam} ${fasta.id}) \\
        | samtools view -bo ${fasta.id}.bam -
    """
}
