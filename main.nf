#!/usr/bin/env nextflow
/*
========================================================================================
                                M e t h y l S c o r e
========================================================================================
 Nextflow implementation of Computomics' MethylScore Pipeline
 #### Homepage / Documentation
 https://github.com/Computomics/MethylScore
 #### Author
 Jörg Hagmann <joerg.hagmann@computomics.com>
 Patrick Hüther <p.huether@lmu.de>
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

log.info """
===================================================================================================================================

███╗   ███╗███████╗████████╗██╗  ██╗██╗   ██╗██╗     ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
████╗ ████║██╔════╝╚══██╔══╝██║  ██║╚██╗ ██╔╝██║     ██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
██╔████╔██║█████╗     ██║   ███████║ ╚████╔╝ ██║     ███████╗██║     ██║   ██║██████╔╝█████╗  
██║╚██╔╝██║██╔══╝     ██║   ██╔══██║  ╚██╔╝  ██║     ╚════██║██║     ██║   ██║██╔══██╗██╔══╝  
██║ ╚═╝ ██║███████╗   ██║   ██║  ██║   ██║   ███████╗███████║╚██████╗╚██████╔╝██║  ██║███████╗
╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝   ╚═╝   ╚══════╝╚══════╝ ╚═════╝ ╚═════╝ ╚═╝  ╚═╝╚══════╝v${workflow.manifest.version}

====================================================================================================================================
Current home              : $HOME
Current user              : $USER
Current path              : $PWD
Script dir                : $projectDir
Working dir               : $workDir
------------------------------------------------------------------------------------------------------------------------------------
PROJECT_FOLDER            : ${params.PROJECT_FOLDER}
REMOVE_INTMED_FILES       : ${params.REMOVE_INTMED_FILES}
------------------------------------------------------------------------------------------------------------------------------------
GENOME                    : ${params.GENOME}
SAMPLE_SHEET              : ${params.SAMPLE_SHEET}
------------------------------------------------------------------------------------------------------------------------------------
BEDGRAPH                  : ${params.BEDGRAPH}
MATRIX                    : ${params.MATRIX}
DO_DEDUP                  : ${(params.BEDGRAPH || params.MATRIX) ? "ignored (BEDGRAPH = ${params.BEDGRAPH})" : params.DO_DEDUP }
STATISTICS                : ${(params.BEDGRAPH || params.MATRIX) ? "ignored (BEDGRAPH = ${params.BEDGRAPH})" : params.STATISTICS }
HUMAN                     : ${params.HUMAN}
IGV OUTPUT                : ${params.IGV}
ROI                       : ${params.ROI}
MR_PARAMS                 : ${params.MR_PARAMS}
METRICS                   : ${params.METRICS}
------------------------------------------------------------------------------------------------------------------------------------
PAIRWISE                  : ${params.PAIRWISE}
DMRS_PER_CONTEXT          : ${params.DMRS_PER_CONTEXT}
DMR_CONTEXTS              : ${params.DMRS_PER_CONTEXT  ? params.DMR_CONTEXTS : 'combined'}
CLUSTER_MIN_METH          : ${params.DMRS_PER_CONTEXT  ? "ignored (DMRS_PER_CONTEXT = ${params.DMRS_PER_CONTEXT})" : params.CLUSTER_MIN_METH}
CLUSTER_MIN_METH_DIFF     : ${params.DMRS_PER_CONTEXT  ? "ignored (DMRS_PER_CONTEXT = ${params.DMRS_PER_CONTEXT})" : params.CLUSTER_MIN_METH_DIFF}
CLUSTER_MIN_METH_CG       : ${!params.DMRS_PER_CONTEXT ? "ignored (DMRS_PER_CONTEXT = ${params.DMRS_PER_CONTEXT})" : params.CLUSTER_MIN_METH_CG}
CLUSTER_MIN_METH_CHG      : ${!params.DMRS_PER_CONTEXT ? "ignored (DMRS_PER_CONTEXT = ${params.DMRS_PER_CONTEXT})" : params.CLUSTER_MIN_METH_CHG}
CLUSTER_MIN_METH_CHH      : ${!params.DMRS_PER_CONTEXT ? "ignored (DMRS_PER_CONTEXT = ${params.DMRS_PER_CONTEXT})" : params.CLUSTER_MIN_METH_CHH}
CLUSTER_MIN_METH_DIFF_CG  : ${!params.DMRS_PER_CONTEXT ? "ignored (DMRS_PER_CONTEXT = ${params.DMRS_PER_CONTEXT})" : params.CLUSTER_MIN_METH_DIFF_CG}
CLUSTER_MIN_METH_DIFF_CHG : ${!params.DMRS_PER_CONTEXT ? "ignored (DMRS_PER_CONTEXT = ${params.DMRS_PER_CONTEXT})" : params.CLUSTER_MIN_METH_DIFF_CHG}
CLUSTER_MIN_METH_DIFF_CHH : ${!params.DMRS_PER_CONTEXT ? "ignored (DMRS_PER_CONTEXT = ${params.DMRS_PER_CONTEXT})" : params.CLUSTER_MIN_METH_DIFF_CHH}
DESERT_SIZE               : ${params.DESERT_SIZE}
DMR_MIN_C                 : ${params.DMR_MIN_C}
DMR_MIN_COV               : ${params.DMR_MIN_COV}
FDR_CUTOFF                : ${params.FDR_CUTOFF}
HDMR_FOLD_CHANGE          : ${params.HDMR_FOLD_CHANGE}
IGNORE_OT                 : ${(params.BEDGRAPH || params.MATRIX) ? "ignored" : params.IGNORE_OT}
IGNORE_OB                 : ${(params.BEDGRAPH || params.MATRIX) ? "ignored" : params.IGNORE_OB}
MERGE_DIST                : ${params.MERGE_DIST}
MIN_COVERAGE              : ${params.MIN_COVERAGE}
MIN_QUAL                  : ${(params.BEDGRAPH || params.MATRIX) ? "ignored" : params.MIN_QUAL }
MR_BATCH_SIZE             : ${params.MR_BATCH_SIZE}
MR_FREQ_CHANGE            : ${params.MR_FREQ_CHANGE}
MR_FREQ_DISTANCE          : ${params.MR_FREQ_DISTANCE}
MR_MIN_C                  : ${params.MR_MIN_C}
SLIDING_WINDOW_SIZE       : ${params.SLIDING_WINDOW_SIZE}
SLIDING_WINDOW_STEP       : ${params.SLIDING_WINDOW_STEP}
TRIM_METHRATE             : ${params.TRIM_METHRATE}
------------------------------------------------------------------------------------------------------------------------------------
Config Profile : ${workflow.profile}
====================================================================================================================================
""".stripIndent()


// validate parameters
ParameterChecks.checkParams(params)

// modules/workflow
if (params.BEDGRAPH) {
    include { BEDGRAPH as CONSENSUS } from './modules/workflow/from_bedgraph'
} else if (params.MATRIX) {
    include { MATRIX as CONSENSUS   } from './modules/workflow/from_matrix'
} else {
    include { BAM as CONSENSUS      } from './modules/workflow/from_alignments'
}

include { SAMPLESHEET               } from './modules/workflow/get_sheets'
include { MRS                       } from './modules/workflow/get_MRs'
include { DMRS                      } from './modules/workflow/get_DMRs'

// modules/process
include { IGV } from './modules/process/generate_igv'

workflow {

    CONSENSUS()

    SAMPLESHEET(CONSENSUS.out.matrixWG)

    def matrix = params.MR_PARAMS ? CONSENSUS.out.matrixCHROM : CONSENSUS.out.matrixWG

    MRS(
        SAMPLESHEET.out.indexedSamples,
        matrix,
        SAMPLESHEET.out.sheet
    )

    if (params.IGV) { IGV(matrixWG, MRS.out.mrs.collect()) }

    DMRS(
        MRS.out.chunks,
        matrix
    )
}
