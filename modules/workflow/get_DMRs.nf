include { CALL_DMRS   } from '../process/call_DMRs'
include { MERGE_DMRS  } from '../process/merge_DMRs'
include { INDEX } from '../process/index_genome_matrix'

workflow DMRS {
    take:
    chunks
    matrix

    main:
    def contexts = params.DMRS_PER_CONTEXT ? Channel.fromList(params.DMR_CONTEXTS.tokenize(',')) : Channel.of('combined')

    INDEX(matrix)

    CALL_DMRS(
        chunks.combine(INDEX.out.index, by:0).transpose(by:2),
        contexts
    )

    MERGE_DMRS(
        CALL_DMRS.out.segments.groupTuple(by:[1,2])
    )

    emit:
    dmrs = MERGE_DMRS.out.dmrs
}