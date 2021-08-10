process IGV {
    tag "${bed}"
    label "resource_low"
    publishDir "${params.PROJECT_FOLDER}/igv", mode: 'copy'

    input:
    tuple val(chromID), path(matrixWG)
    path(bed)

    output:
    path('methinfo.igv'), emit: igv

    script:
    """
    sort -m -k1,1 -k2,2g -k3,3g ${bed} > MRs.merged.bed
    python matrix2igv.py  \\
        -i ${matrixWG} \\
        -m MRs.merged.bed \\
        -o methinfo.igv
    """
}