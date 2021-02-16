include { GET_SAMPLES     } from './get_samples'
include { MATRIX          } from './get_genome_matrix'
include { SORT_BAM        } from '../process/sort_bam_samtools'
include { MERGE_BAM       } from '../process/MergeSamFiles_picard'
include { DEDUPLICATE     } from '../process/MarkDuplicates_picard'
include { READ_STATISTICS } from '../process/get_read_statistics'
include { SPLIT_BAM       } from '../process/split_bam_by_chromosome'
include { METHYLDACKEL    } from '../process/methyldackel'

workflow BAM {
    main:

    def roi_file = params.ROI ? Channel.fromPath(params.ROI, checkIfExists: true).collect() : file('null')

    GET_SAMPLES()

    SORT_BAM(GET_SAMPLES.out.samples)

    SORT_BAM.out.bam
        .groupTuple(by:0)
        .branch { sampleID, bam ->
            multi_rep: bam.toList().size() > 1
            single_rep: true
        }
        .set { samples }

    MERGE_BAM(samples.multi_rep)

    def alignments = params.DO_DEDUP ? samples.single_rep | mix(MERGE_BAM.out.bam) | DEDUPLICATE : samples.single_rep | mix(MERGE_BAM.out.bam)

    if ( params.STATISTICS ) { READ_STATISTICS(alignments, roi_file) }

    SPLIT_BAM(alignments, GET_SAMPLES.out.fasta)

    METHYLDACKEL(SPLIT_BAM.out.bam)

    MATRIX(
        METHYLDACKEL.out.consensus,
        alignments
    )

    emit:
    matrixWG = MATRIX.out.matrixWG
    matrixCHROM = MATRIX.out.matrixCHROM
}