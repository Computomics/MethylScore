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


// modules/process
include { CALL_MRS   } from './modules/process/call_MRs'
include { SPLIT_MRS  } from './modules/process/split_MRs'
include { CALL_DMRS  } from './modules/process/call_DMRs'
include { MERGE_DMRS } from './modules/process/merge_DMRs'

// modules/workflow
if (params.BEDGRAPH) {
    include { BEDGRAPH as CONSENSUS } from './modules/workflow/from_bedgraph'
} else if (params.MATRIX) {
    include { MATRIX as CONSENSUS   } from './modules/workflow/from_matrix'
} else {
    include { BAM as CONSENSUS      } from './modules/workflow/from_alignments'
}

workflow {

    contexts = params.DMRS_PER_CONTEXT ? Channel.fromList(params.DMR_CONTEXTS.tokenize(',')) : Channel.of('combined')
    hmm_params_file = params.MR_PARAMS ? Channel.fromPath(params.MR_PARAMS, checkIfExists: true).collect() : file('null')

    CONSENSUS()

    if (params.MR_PARAMS) {
        CONSENSUS.out.matrixCHROM.set{ matrix }
    } else {
        CONSENSUS.out.matrixWG.set{ matrix }
    }

    CALL_MRS(
        CONSENSUS.out.indexedSamples,
        matrix,
        hmm_params_file
    )

    if (params.IGV) { MATRIX_TO_IGV(CONSENSUS.out.matrixWG, CALL_MRS.out.bed.collect{ it[1] }) }

    SPLIT_MRS(
        CALL_MRS.out.bed.collectFile(cache:true){ sample, bed -> ["${sample}.MRs.bed", bed] }.collect(),
        CONSENSUS.out.mrsheet
    )

    CALL_DMRS(
        CONSENSUS.out.index.collect(),
        CONSENSUS.out.dmrsheet.collect(),
        SPLIT_MRS.out.chunks.transpose(),
        contexts
    )

    MERGE_DMRS(
        CALL_DMRS.out.segments.collectFile(cache:true){ comp, context, segment -> ["${comp}.${context}.dif", segment] },
        CONSENSUS.out.dmrsheet
    )
}
