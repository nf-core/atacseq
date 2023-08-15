include { SAMTOOLS_SORT as SAMTOOLS_SORT_SE } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_PE } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX          } from '../../modules/nf-core/samtools/index/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../nf-core/bam_sort_stats_samtools/main'
include { BAM_STATS_SAMTOOLS      } from '../nf-core/bam_stats_samtools/main'

workflow BAM_NOFILTER_BAMTOOLS {
    take:
    ch_bam                  // channel: [ val(meta), [ bam ] ]
    ch_fasta                // channel: [ fasta ]

    main:

    ch_versions = Channel.empty()

    ch_bam
        .branch {
            meta, bam ->
                single_end: meta.single_end
                    return [ meta, bam ]
                paired_end: !meta.single_end
                    return [ meta, bam ]
        }
        .set { ch_bam_b }

    //
    // Index SE BAM file
    //
    SAMTOOLS_INDEX {
        ch_bam_b.single_end
    }
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    //
    // Name sort SE BAM
    //
    SAMTOOLS_SORT_SE (
        ch_bam_b.single_end
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_SE.out.versions.first())

    //
    // Run samtools stats, flagstat and idxstats on SE BAM
    //
    BAM_STATS_SAMTOOLS (
        ch_bam_b.single_end.join(SAMTOOLS_INDEX.out.bai),
        ch_fasta
    )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions.first())

    //
    // Name sort PE BAM
    //
    SAMTOOLS_SORT_PE (
        ch_bam_b.paired_end
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_PE.out.versions.first())

    //
    // Sort, index PE BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS (
        SAMTOOLS_SORT_PE.out.bam,
        ch_fasta
    )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions.first())

    emit:
    name_bam = SAMTOOLS_SORT_PE.out.bam.mix(SAMTOOLS_SORT_SE.out.bam)                            // channel: [ val(meta), [ bam ] ]
    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam.mix(ch_bam_b.single_end)                  // channel: [ val(meta), [ bam ] ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats.mix(BAM_STATS_SAMTOOLS.out.stats)       // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat.mix(BAM_STATS_SAMTOOLS.out.flagstat) // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats.mix(BAM_STATS_SAMTOOLS.out.idxstats) // channel: [ val(meta), [ idxstats ] ]
    versions = ch_versions                                                               // channel: [ versions.yml ]
}
