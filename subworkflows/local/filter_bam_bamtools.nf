include { BAMTOOLS_FILTER         } from '../../modules/local/bamtools_filter'
include { BAM_REMOVE_ORPHANS      } from '../../modules/local/bam_remove_orphans'
include { BAM_SORT_STATS_SAMTOOLS } from '../nf-core/bam_sort_stats_samtools/main'

workflow FILTER_BAM_BAMTOOLS {
    take:
    ch_bam_bai                   // channel: [ val(meta), [ bam ], [bai] ]
    ch_bed                       // channel: [ bed ]
    ch_fasta                     // channel: [ fasta ]
    ch_bamtools_filter_se_config // channel: [ config_file ]
    ch_bamtools_filter_pe_config // channel: [ config_file ]

    main:

    ch_versions = Channel.empty()

    //
    // Filter BAM file with BAMTools
    //
    BAMTOOLS_FILTER (
        ch_bam_bai,
        ch_bed,
        ch_bamtools_filter_se_config,
        ch_bamtools_filter_pe_config
    )
    ch_versions = ch_versions.mix(BAMTOOLS_FILTER.out.versions.first())

    //
    // Remove orphan reads from BAM file
    //
    BAM_REMOVE_ORPHANS (
        BAMTOOLS_FILTER.out.bam
    )
    ch_versions = ch_versions.mix(BAM_REMOVE_ORPHANS.out.versions.first())

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS (
        BAM_REMOVE_ORPHANS.out.bam,
        ch_fasta
    )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions.first())

    emit:
    name_bam = BAM_REMOVE_ORPHANS.out.bam           // channel: [ val(meta), [ bam ] ]

    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai      = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                          // channel: [ versions.yml ]
}
