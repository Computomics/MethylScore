include { GET_SAMPLES } from './sub/get_samples'
include { SORT        } from '../process/sort_bam_samtools'
include { MERGE       } from '../process/MergeSamFiles_picard'
include { DEDUP       } from '../process/MarkDuplicates_picard'
include { STATISTICS  } from '../process/get_read_statistics'
include { SPLIT       } from '../process/split_bam_by_chromosome'
include { MATRIX      } from './sub/get_genome_matrix'

include { METHYLDACKEL } from '../process/methyldackel'

workflow BAM {
    main:

    roi_file = params.ROI ? Channel.fromPath(params.ROI, checkIfExists: true).collect() : file('null')

    GET_SAMPLES()

    SORT(GET_SAMPLES.out.samples)

    SORT.out.bam
        .groupTuple(by:0)
        .branch { sampleID, bam ->
            multi_rep: bam.toList().size() > 1
            single_rep: true
        }
        .set { samples }

    samples.multi_rep | MERGE

    alignments = samples.single_rep.mix(MERGE.out.bam)

    if ( params.DO_DEDUP ) { alignments = DEDUP(alignments).bam }

    if ( params.STATISTICS ) { STATISTICS(alignments, roi_file) }

    SPLIT(alignments, GET_SAMPLES.out.fasta)

    METHYLDACKEL(SPLIT.out.bam.join(SPLIT.out.fasta))

    MATRIX(
        METHYLDACKEL.out.consensus,
        SPLIT.out.fasta,
        alignments
    )

    emit:
    matrixWG = MATRIX.out.matrixWG
    indexedSamples = MATRIX.out.indexedSamples
    mrsheet = MATRIX.out.mrsheet
    dmrsheet = MATRIX.out.dmrsheet
    index = MATRIX.out.index
}