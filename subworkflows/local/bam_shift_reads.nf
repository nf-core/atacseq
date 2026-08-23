include { SAMTOOLS_SORT            } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX           } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT        } from '../../modules/nf-core/samtools/flagstat/main'
include { DEEPTOOLS_ALIGNMENTSIEVE } from '../../modules/nf-core/deeptools/alignmentsieve'

workflow BAM_SHIFT_READS {
    take:
    ch_bam_index // channel: [ val(meta), [ bam ], [ bai/csi ] ]
    ch_fasta_fai // channel: [ val(meta), path(fasta), path(fai) ]

    main:

    //
    // Shift reads
    //
    DEEPTOOLS_ALIGNMENTSIEVE (
        ch_bam_index
    )

    //
    // Sort reads
    //
    SAMTOOLS_SORT (
        DEEPTOOLS_ALIGNMENTSIEVE.out.bam,
        ch_fasta_fai,
        ''
    )

    //
    // Index reads
    //
    SAMTOOLS_INDEX (
        SAMTOOLS_SORT.out.bam
    )

    //
    // Run samtools flagstat
    //
    SAMTOOLS_FLAGSTAT (
        SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.index, by: [0])
    )

    emit:
    bam      = SAMTOOLS_SORT.out.bam                // channel: [ val(meta), [ bam ] ]
    index    = SAMTOOLS_INDEX.out.index             // channel: [ val(meta), [ bai/csi ] ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat       // channel: [ val(meta), [ flagstat ] ]
}
