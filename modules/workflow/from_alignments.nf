include { GET_SAMPLES } from './get_samples'
include { MATRIX      } from './get_genome_matrix'
include { SORT        } from '../process/sort_bam_samtools'
include { MERGE       } from '../process/MergeSamFiles_picard'
include { DEDUP       } from '../process/MarkDuplicates_picard'
include { STATISTICS  } from '../process/get_read_statistics'
include { SPLIT       } from '../process/split_bam_by_chromosome'

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

    METHYLDACKEL(SPLIT.out.bam)

    MATRIX(
        METHYLDACKEL.out.consensus,
        alignments
    )

    emit:
    matrixWG = MATRIX.out.matrixWG
    matrixCHROM = MATRIX.out.matrixCHROM
}