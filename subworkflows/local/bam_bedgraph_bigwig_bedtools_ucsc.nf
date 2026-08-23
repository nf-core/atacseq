
//
// Convert BAM to normalised bigWig via bedGraph using BEDTools and UCSC
//

include { BEDTOOLS_GENOMECOV    } from '../../modules/local/bedtools_genomecov'
include { UCSC_BEDGRAPHTOBIGWIG } from '../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow BAM_BEDGRAPH_BIGWIG_BEDTOOLS_UCSC {
    take:
    ch_bam_flagstat // channel: [ val(meta), [bam], [flagstat] ]
    ch_chrom_sizes  // channel: [ bed ]

    main:


    //
    // Create bedGraph coverage track
    //
    BEDTOOLS_GENOMECOV (
        ch_bam_flagstat
    )

    //
    // Create bigWig coverage tracks
    //
    UCSC_BEDGRAPHTOBIGWIG (
        BEDTOOLS_GENOMECOV.out.bedgraph,
        ch_chrom_sizes
    )

    emit:
    bedgraph     = BEDTOOLS_GENOMECOV.out.bedgraph      // channel: [ val(meta), [ bedgraph ] ]
    scale_factor = BEDTOOLS_GENOMECOV.out.scale_factor  // channel: [ val(meta), [ txt ] ]

    bigwig       = UCSC_BEDGRAPHTOBIGWIG.out.bigwig     // channel: [ val(meta), [ bigwig ] ]

}
