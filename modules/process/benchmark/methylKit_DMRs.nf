process RUN_METHYLKIT {
    tag "${fasta.id}"
    publishDir "${params.PROJECT_FOLDER}/methylKit/dmrs", mode: 'copy'
    container 'quay.io/biocontainers/bioconductor-methylkit:1.16.1--r40h399db7b_0'

    input:
    tuple val(fasta), val(group), val(sampleID), val(context), path(inputfile)

    output:
    path("*.DMRs.bed"), emit: methylkit_dmrs

    script:
    """
    #!/usr/bin/env Rscript

    library(methylKit)
    
    methobj <- methRead(list("${inputfile.join('","')}"),
                        sample.id=list("${sampleID.join('","')}"),
                        assembly="",
                        treatment=c(${sampleID.collect{ it.startsWith(group.unique()[0]) ? 1 : 0 }.join(',')}),
                        context="${context}",
                        mincov=${params.DMR_MIN_COV})

    tiles = tileMethylCounts(methobj,win.size=200,step.size=200,cov.bases=${params.DMR_MIN_COV})

    meth <- unite(tiles,destrand=FALSE)
    meth.diff <- calculateDiffMeth(meth,mc.cores=${task.cpus})
    dmrs <- data.frame(getMethylDiff(meth.diff,difference=${params.CLUSTER_MIN_METH_DIFF},qvalue=${params.FDR_CUTOFF}))
    
    write.table(dmrs, file="methylKit.${context}.DMRs.bed", quote=FALSE, sep='\\t', row.names = F) 
    """
}
