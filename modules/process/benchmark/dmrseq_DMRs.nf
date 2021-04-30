process RUN_DMRSEQ {
    tag "${context}"
    publishDir "${params.PROJECT_FOLDER}/dmrseq/dmrs", mode: 'copy'
    container 'quay.io/biocontainers/bioconductor-dmrseq:1.10.0--r40_0'

    input:
    tuple val(group), val(sampleID), val(context), path(inputfile), val(fasta)

    output:
    path("*.DMRs.bed"), emit: dmrseq_dmrs

    script:
    """
    #!/usr/bin/env Rscript

    library(BiocParallel)
    library(dmrseq)
    library(dplyr)
    library(readr)

    register(MulticoreParam(${task.cpus}))

    data <- read.bismark(c("${inputfile.join('","')}"))

    pData(data)\$group <- c("${group.join('","')}")
    pData(data)\$rep <- c("${sampleID.join('","')}")

    dmrs <- dmrseq(data, cutoff=${params.FDR_CUTOFF}, testCovariate="group")
    
    dmrs %>% as_tibble() %>% arrange(start) %>% write_tsv("dmrseq.${context}.DMRs.bed")
    """
}