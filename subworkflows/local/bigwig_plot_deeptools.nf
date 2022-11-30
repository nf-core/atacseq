
//
// Create coverage matrices, plot coverage profiles and heatmap with deepTools
//

include { DEEPTOOLS_PLOTPROFILE } from '../../modules/nf-core/deeptools/plotprofile/main'
include { DEEPTOOLS_PLOTHEATMAP } from '../../modules/nf-core/deeptools/plotheatmap/main'
include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_SCALE_REGIONS   } from '../../modules/nf-core/deeptools/computematrix/main'
include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_REFERENCE_POINT } from '../../modules/nf-core/deeptools/computematrix/main'

workflow BIGWIG_PLOT_DEEPTOOLS {
    take:
    ch_bigwig   // channel: [ val(meta), bigwig ]
    ch_gene_bed // channel: [ bed ]
    ch_tss_bed  // channel: [ bed ]
    
    main:

    ch_versions = Channel.empty()

    //
    // deepTools matrix generation for plotting over full transcript length
    //
    DEEPTOOLS_COMPUTEMATRIX_SCALE_REGIONS (
        ch_bigwig,
        ch_gene_bed
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_COMPUTEMATRIX_SCALE_REGIONS.out.versions.first())

    //
    // deepTools matrix generation for plotting at TSS point
    //
    DEEPTOOLS_COMPUTEMATRIX_REFERENCE_POINT (
        ch_bigwig,
        ch_tss_bed
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_COMPUTEMATRIX_REFERENCE_POINT.out.versions.first())

    //
    // deepTools profile plots
    //
    DEEPTOOLS_PLOTPROFILE (
        DEEPTOOLS_COMPUTEMATRIX_SCALE_REGIONS.out.matrix
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_PLOTPROFILE.out.versions.first())

    //
    // deepTools heatmaps
    //
    DEEPTOOLS_PLOTHEATMAP (
        DEEPTOOLS_COMPUTEMATRIX_REFERENCE_POINT.out.matrix
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_PLOTHEATMAP.out.versions.first())

    emit:
    scale_regions_matrix   = DEEPTOOLS_COMPUTEMATRIX_SCALE_REGIONS.out.matrix   // channel: [ val(meta), [ matrix ] ]
    scale_regions_table    = DEEPTOOLS_COMPUTEMATRIX_SCALE_REGIONS.out.table    // channel: [ val(meta), [ table ] ]

    reference_point_matrix = DEEPTOOLS_COMPUTEMATRIX_REFERENCE_POINT.out.matrix // channel: [ val(meta), [ matrix ] ]
    reference_point_table  = DEEPTOOLS_COMPUTEMATRIX_REFERENCE_POINT.out.table  // channel: [ val(meta), [ table ] ]

    plotprofile_pdf        = DEEPTOOLS_PLOTPROFILE.out.pdf                      // channel: [ val(meta), [ pdf ] ]
    plotprofile_table      = DEEPTOOLS_PLOTPROFILE.out.table                    // channel: [ val(meta), [ table ] ]

    plotheatmap_pdf        = DEEPTOOLS_PLOTHEATMAP.out.pdf                      // channel: [ val(meta), [ pdf ] ]
    plotheatmap_table      = DEEPTOOLS_PLOTHEATMAP.out.table                    // channel: [ val(meta), [ table ] ]

    versions               = ch_versions                                        // channel: [ versions.yml ]
}
