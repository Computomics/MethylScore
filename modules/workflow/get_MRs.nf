include { CALL_MRS      } from '../process/call_MRs'
include { MR_STATISTICS } from '../process/get_MR_statistics'
include { SPLIT_MRS     } from '../process/split_MRs'

workflow MRS {
    take:
    samples
    matrix
    samplesheet

    main:
    def hmm_params_file = params.MR_PARAMS ? Channel.fromPath(params.MR_PARAMS, checkIfExists: true).collect() : file('null')

    CALL_MRS(
        samples.combine(matrix),
        hmm_params_file
    )

    MR_STATISTICS( CALL_MRS.out.bed.collectFile(cache:true, storeDir:"${params.PROJECT_FOLDER}/04MRs", sort: { it[0] }){ chrom, sample, bed -> ["${sample}.MRs.bed", bed] } )
 
    SPLIT_MRS(
        CALL_MRS.out.bed.groupTuple(by:0),
        samplesheet
    )

    emit:
    chunks = SPLIT_MRS.out.chunks
}