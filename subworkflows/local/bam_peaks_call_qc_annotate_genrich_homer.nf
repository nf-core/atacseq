//
// Call peaks with Genrich, annotate with HOMER and perform downstream QC
//

include { GENRICH           } from '../../modules/nf-core/genrich/main'
include { HOMER_ANNOTATEPEAKS      } from '../../modules/nf-core/homer/annotatepeaks/main'

include { FRIP_SCORE               } from '../../modules/local/frip_score'
include { MULTIQC_CUSTOM_PEAKS     } from '../../modules/local/multiqc_custom_peaks'

include { PLOT_PEAKS_QC as PLOT_GENRICH_QC      } from '../../modules/local/plot_peaks_qc'

include { PLOT_HOMER_ANNOTATEPEAKS } from '../../modules/local/plot_homer_annotatepeaks'

workflow BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER {
    take:
    ch_bam                            // channel: [ val(meta), [ treatment_bam ], [ control_bam ] ]
    ch_fasta                          // channel: [ fasta ]
    ch_gtf                            // channel: [ gtf ]
    genome                            // string: genome name
    ch_blacklist_regions              // channel: [ bed ]
    annotate_peaks_suffix             // string: suffix for input HOMER annotate peaks files to be trimmed off
    ch_peak_count_header_multiqc      // channel: [ header_file ]
    ch_frip_score_multiqc             // channel: [ header_file ]
    ch_peak_annotation_header_multiqc // channel: [ header_file ]
// add broad_peak options
    is_narrow_peak                    // boolean: true/false
    skip_peak_annotation              // boolean: true/false
    skip_peak_qc                      // boolean: true/false

    save_pvalues                      // boolean: true/false
    save_pileup                       // boolean: true/false
    save_bed                          // boolean: true/false
    save_duplicates                   // boolean: true/false

    main:

    ch_versions = Channel.empty()

    //
    // Call peaks with Genrich
    //
    GENRICH (
        ch_bam,
        ch_blacklist_regions,
        save_pvalues,
        save_pileup,
        save_bed,
        save_duplicates

    )
    ch_versions = ch_versions.mix(GENRICH.out.versions.first())

    //
    // Filter out samples with 0 Genrich peaks called
    //
    GENRICH
        .out
        .peak
        .filter {
            meta, peaks ->
                peaks.size() > 0
        }
        .set { ch_genrich_peaks }


    // Create channels: [ meta, ip_bam, peaks ]

    ch_bam
        .join(ch_genrich_peaks, by: [0])
        .map {
            meta, treatment_bam, control_bam, peaks ->
                [ meta, treatment_bam, peaks ]
        }
        .set { ch_bam_peaks }

    // split ch_bam_peaks by treatment_bam

    ch_bam_peaks
        .transpose()
        .map {
            meta, treatment_bam, peaks ->
                def meta_clone = meta.clone()
                meta_clone.id = treatment_bam.getSimpleName()
                [ meta_clone, treatment_bam, peaks ]
        }
        .set { ch_bam_peaks_treatment_bam }

    //
    // Calculate FRiP score
    //
    FRIP_SCORE (
        ch_bam_peaks_treatment_bam
    )
    ch_versions = ch_versions.mix(FRIP_SCORE.out.versions.first())

    // Create channels: [ meta, peaks, frip ]
    ch_bam_peaks_treatment_bam
        .join(FRIP_SCORE.out.txt, by: [0])
        .map {
            meta, ip_bam, peaks, frip ->
                [ meta, peaks, frip ]
        }
        .set { ch_bam_peak_frip }

    //
    // FRiP score custom content for MultiQC
    //
    MULTIQC_CUSTOM_PEAKS (
        ch_bam_peak_frip,
        ch_peak_count_header_multiqc,
        ch_frip_score_multiqc
    )
    ch_versions = ch_versions.mix(MULTIQC_CUSTOM_PEAKS.out.versions.first())

    ch_homer_annotatepeaks              = Channel.empty()
    ch_homer_det_annotatepeaks          = Channel.empty()
    ch_plot_genrich_qc_txt              = Channel.empty()
    ch_plot_genrich_qc_pdf              = Channel.empty()
    ch_plot_homer_annotatepeaks_txt     = Channel.empty()
    ch_plot_homer_annotatepeaks_pdf     = Channel.empty()
    ch_plot_homer_annotatepeaks_tsv     = Channel.empty()
    if (!skip_peak_annotation) {
        //
        // Annotate peaks with HOMER
        //
        HOMER_ANNOTATEPEAKS (
            ch_genrich_peaks,
            ch_fasta,
            ch_gtf
        )
        ch_homer_annotatepeaks = HOMER_ANNOTATEPEAKS.out.txt
        ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS.out.versions.first())

        if (!skip_peak_qc) {
            //
            // Genrich QC plots with R
            //
            PLOT_GENRICH_QC (
                ch_genrich_peaks.collect{it[1]},
                is_narrow_peak
            )
            ch_plot_genrich_qc_txt = PLOT_GENRICH_QC.out.txt
            ch_plot_genrich_qc_pdf = PLOT_GENRICH_QC.out.pdf
            ch_versions = ch_versions.mix(PLOT_GENRICH_QC.out.versions)

            //
            // Peak annotation QC plots with R
            //
            PLOT_HOMER_ANNOTATEPEAKS (
                HOMER_ANNOTATEPEAKS.out.txt.collect{it[1]},
                ch_peak_annotation_header_multiqc,
                annotate_peaks_suffix
            )
            ch_plot_homer_annotatepeaks_txt = PLOT_HOMER_ANNOTATEPEAKS.out.txt
            ch_plot_homer_annotatepeaks_pdf = PLOT_HOMER_ANNOTATEPEAKS.out.pdf
            ch_plot_homer_annotatepeaks_tsv = PLOT_HOMER_ANNOTATEPEAKS.out.tsv
            ch_versions = ch_versions.mix(PLOT_HOMER_ANNOTATEPEAKS.out.versions)
        }
    }

    emit:
    peaks                        = ch_genrich_peaks                   // channel: [ val(meta), [ peaks ] ]
    bed_intervals                = GENRICH.out.bed_intervals           // channel: [ val(meta), [ bed ] ]
    bedgraph_pileup              = GENRICH.out.bedgraph_pileup           // channel: [ val(meta), [ bedgraph ] ]
    bedgraph_pvalues             = GENRICH.out.bedgraph_pvalues           // channel: [ val(meta), [ bedgraph ] ]
    duplicates                   = GENRICH.out.duplicates           // channel: [ val(meta), [ bedgraph ] ]

    frip_txt                     = FRIP_SCORE.out.txt               // channel: [ val(meta), [ txt ] ]

    frip_multiqc                 = MULTIQC_CUSTOM_PEAKS.out.frip    // channel: [ val(meta), [ frip ] ]
    peak_count_multiqc           = MULTIQC_CUSTOM_PEAKS.out.count   // channel: [ val(meta), [ counts ] ]

    homer_annotatepeaks          = ch_homer_annotatepeaks           // channel: [ val(meta), [ txt ] ]
    plot_genrich_qc_txt          = ch_plot_genrich_qc_txt             // channel: [ txt ]
    plot_genrich_qc_pdf          = ch_plot_genrich_qc_pdf             // channel: [ pdf ]

    plot_homer_annotatepeaks_txt = ch_plot_homer_annotatepeaks_txt  // channel: [ txt ]
    plot_homer_annotatepeaks_pdf = ch_plot_homer_annotatepeaks_pdf  // channel: [ pdf ]
    plot_homer_annotatepeaks_tsv = ch_plot_homer_annotatepeaks_tsv  // channel: [ tsv ]

    versions                     = ch_versions                      // channel: [ versions.yml ]
}
