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

    CALL_DMRS.out.segments
        .collectFile(cache:true){ comp, context, segment -> ["${comp}__${context}.dif", segment] }
        .map { segments -> def keys = segments.baseName.minus('.dif').split('__'); [keys[0], keys[1], segments] }
        .join(CALL_DMRS.out.sheet, by:[0,1])
        .set { segments }

    MERGE_DMRS(segments)

    emit:
    dmrs = MERGE_DMRS.out.dmrs
}