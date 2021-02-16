include { CALL_DMRS   } from '../process/call_DMRs'
include { MERGE_DMRS  } from '../process/merge_DMRs'
include { BUILD_INDEX } from '../process/index_genome_matrix'

workflow DMRS {
    take:
    chunks
    matrix
    samplesheet

    main:
    def contexts = params.DMRS_PER_CONTEXT ? Channel.fromList(params.DMR_CONTEXTS.tokenize(',')) : Channel.of('combined')

    BUILD_INDEX(matrix)

    CALL_DMRS(
        chunks.combine(BUILD_INDEX.out.index, by:0).transpose(by:2),
        contexts
    )

    MERGE_DMRS(
        CALL_DMRS.out.segments.collectFile(cache:true){ comp, context, segment -> ["${comp}.${context}.dif", segment] },
        samplesheet
    )

    emit:
    dmrs = MERGE_DMRS.out.dmrs
}