process RUN_DSS {
    tag "${fasta.id}"
    publishDir "${params.PROJECT_FOLDER}/dss/dmrs", mode: 'copy'
    container 'quay.io/biocontainers/bioconductor-dss:2.38.0--r40h037d062_0'

    input:
    tuple val(fasta), val(group), val(sampleID), val(context), path(inputfile)

    output:
    path("*.DMRs.bed"), emit: dss_dmrs

    script:
    """
    #!/usr/bin/env Rscript

    library(DSS)

    multicore <- MulticoreParam(workers=${task.cpus})

    headerline <- c("chr","pos","N","X")
    samples <- c("${sampleID.join('","')}")
    files <- c("${inputfile.join('","')}")

    for (i in 1:length(files)) {
        assign(samples[i],read.table(files[i],header=F,col.names=headerline))
    }

    BSobj <- makeBSseqData(mget(samples), samples)

    dmlTest <- DMLtest(BSobj,
                       group1=c("${sampleID.findAll{ it.startsWith(group.unique()[0]) }.join('","')}"),
                       group2=c("${sampleID.findAll{ it.startsWith(group.unique()[1]) }.join('","')}"),
                       smoothing=TRUE,
                       BPPARAM = multicore)

    dmrs <- callDMR(dmlTest,
                delta=${params.CLUSTER_MIN_METH_DIFF/100},
                minCG=${params.DMR_MIN_C},
                p.threshold=${params.FDR_CUTOFF})
    
    write.table(dmrs, file="DSS.${context}.DMRs.bed", quote=FALSE, sep='\\t', row.names = F) 
    """
}