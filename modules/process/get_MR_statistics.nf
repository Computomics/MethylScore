process MR_STATISTICS {
    tag "${sampleID}"
    label "resource_low"
    publishDir "${params.PROJECT_FOLDER}/04MRs/stats", mode: 'copy'

    input:
    path(bed)

    output:
    path("${sampleID}.MR_stats.tsv"), emit: stats

    shell:
    sampleID = bed.name.minus('.MRs.bed') 
    '''
    echo -e "#sampleID\tMRcount\tgenomespace\tavg_length" > !{sampleID}.MR_stats.tsv
    echo -e "!{sampleID}\t$(cat !{bed}| wc -l)\t$(awk -v OFS='\t' '{sum+=$3-$2+1}END{print sum, sprintf("%.0f", sum/NR)}' !{bed})" >> !{sampleID}.MR_stats.tsv
    '''
}