include { SAMTOOLS_SORT           } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX          } from '../../modules/nf-core/samtools/index/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../nf-core/bam_sort_stats_samtools/main'
include { BAM_STATS_SAMTOOLS      } from '../nf-core/bam_stats_samtools/main'

include { BAMTOOLS_FILTER         } from '../../modules/local/bamtools_filter'
include { BAM_REMOVE_ORPHANS      } from '../../modules/local/bam_remove_orphans'

workflow BAM_FILTER_BAMTOOLS {
    take:
    ch_bam_index                 // channel: [ val(meta), [ bam ], [ bai/csi ] ]
    ch_bed                       // channel: [ bed ]
    ch_fasta_fai                 // channel: [ val(meta), path(fasta), path(fai) ]
    ch_bamtools_filter_se_config // channel: [ config_file ]
    ch_bamtools_filter_pe_config // channel: [ config_file ]

    main:


    //
    // Filter BAM file with BAMTools
    //
    BAMTOOLS_FILTER (
        ch_bam_index,
        ch_bed,
        ch_bamtools_filter_se_config,
        ch_bamtools_filter_pe_config
    )

    BAMTOOLS_FILTER
        .out
        .bam
        .branch {
            meta, bam ->
                single_end: meta.single_end
                    return [ meta, bam ]
                paired_end: !meta.single_end
                    return [ meta, bam ]
        }
        .set { ch_bam }

    //
    // Index SE BAM file
    //
    SAMTOOLS_INDEX {
        ch_bam.single_end
    }

    ch_index = SAMTOOLS_INDEX.out.index

    //
    // Run samtools stats, flagstat and idxstats on SE BAM
    //
    BAM_STATS_SAMTOOLS (
        ch_bam.single_end.join(ch_index),
        ch_fasta_fai
    )

    //
    // Name sort PE BAM before filtering with pysam
    //
    SAMTOOLS_SORT (
        ch_bam.paired_end,
        ch_fasta_fai,
        ''
    )

    //
    // Remove orphan reads from PE BAM file
    //
    BAM_REMOVE_ORPHANS (
        SAMTOOLS_SORT.out.bam
    )

    //
    // Sort, index PE BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS (
        BAM_REMOVE_ORPHANS.out.bam,
        ch_fasta_fai
    )

    emit:
    name_bam = SAMTOOLS_SORT.out.bam                                                     // channel: [ val(meta), [ bam ] ]
    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam.mix(ch_bam.single_end)                    // channel: [ val(meta), [ bam ] ]
    index    = BAM_SORT_STATS_SAMTOOLS.out.index.mix(ch_index)                           // channel: [ val(meta), [ bai/csi ] ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats.mix(BAM_STATS_SAMTOOLS.out.stats)       // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat.mix(BAM_STATS_SAMTOOLS.out.flagstat) // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats.mix(BAM_STATS_SAMTOOLS.out.idxstats) // channel: [ val(meta), [ idxstats ] ]
}
