//
// Call consensus peaks with BEDTools and custom scripts, annotate with HOMER, quantify with featureCounts and QC with DESeq2
//

include { HOMER_ANNOTATEPEAKS    } from '../../modules/nf-core/homer/annotatepeaks/main'
include { SUBREAD_FEATURECOUNTS  } from '../../modules/nf-core/subread/featurecounts/main'

include { MACS3_CONSENSUS        } from '../../modules/local/macs3_consensus'
include { FEATURECOUNTS_MERGE    } from '../../modules/local/featurecounts_merge'
include { DESEQ2_QC              } from '../../modules/local/deseq2_qc'

workflow BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2 {
    take:
    ch_peaks                            // channel: [ val(meta), [ peaks ] ]
    ch_bams                             // channel: [ val(meta), [ bams ] ]
    ch_fasta                            // channel: [ fasta ]
    ch_gtf                              // channel: [ gtf ]
    ch_deseq2_pca_header_multiqc        // channel: [ header_file ]
    ch_deseq2_clustering_header_multiqc // channel: [ header_file ]
    is_narrow_peak                      // boolean: true/false
    skip_peak_annotation                // boolean: true/false
    skip_deseq2_qc                      // boolean: true/false

    main:


    // Create channels: [ meta , [ peaks ] ]
    // where meta = [ id : consensus_peaks ]
    ch_peaks
        .collect { item -> item[1] }
        .filter { item -> item.size() > 1 }
        .map {
            peaks ->
                [ [ id: 'consensus_peaks' ], peaks ]
        }
        .set { ch_consensus_peaks }

    //
    // Generate consensus peaks across samples
    //
    MACS3_CONSENSUS (
        ch_consensus_peaks,
        is_narrow_peak
    )

    //
    // Annotate consensus peaks
    //
    ch_homer_annotatepeaks = channel.empty()
    if (!skip_peak_annotation) {
        HOMER_ANNOTATEPEAKS (
            MACS3_CONSENSUS.out.bed,
            ch_fasta,
            ch_gtf
        )
        ch_homer_annotatepeaks = HOMER_ANNOTATEPEAKS.out.txt
    }

    //
    // Quantify peaks across samples with featureCounts.
    //
    // featureCounts (subread >= 2.1.0) applies paired-end mode (-p) to a whole
    // invocation and aborts when that invocation mixes single-end and paired-end
    // BAMs. The consensus BAMs can span both library types, so split them by
    // endedness, count each homogeneous batch with the correct pairing flag
    // (derived from meta.single_end inside SUBREAD_FEATURECOUNTS), then merge the
    // per-batch matrices back into one consensus table for DESeq2 and MultiQC.
    // The join with ch_peaks keeps only samples that contributed peaks; combining
    // with MACS3_CONSENSUS.out.saf also gates counting on a consensus existing
    // (>= 2 samples), matching the previous behaviour.
    //
    ch_consensus_saf = MACS3_CONSENSUS.out.saf.map { _meta, saf -> saf }

    // The merged-library caller joins in a control-BAM column
    // ([ meta, bam, control ] -> [ meta, bam, control, peak ]) while the
    // merged-replicate caller does not ([ meta, bams ] -> [ meta, bams, peak ]),
    // so the joined tuple arity differs between the two instantiations of this
    // subworkflow. Index positionally (meta = item[0], bam = item[1]) to stay
    // tolerant of both shapes, as the pre-split implementation did.
    ch_bams
        .join(ch_peaks)
        .branch { item ->
            single_end: item[0].single_end
            paired_end: !item[0].single_end
        }
        .set { ch_consensus_bams }

    // Each batch is assembled from an unordered channel collect, so sort by
    // filename: the BAM order sets the featureCounts column order, and an
    // unsorted list makes the count matrix (and its snapshot md5) vary between
    // runs and hosts.
    ch_se_batch = ch_consensus_bams.single_end
        .map { item -> item[1] }
        .collect()
        .filter { bams -> bams }
        .map { bams -> [ [ id: 'consensus_peaks', single_end: true ], bams.toSorted { it.name } ] }

    ch_pe_batch = ch_consensus_bams.paired_end
        .map { item -> item[1] }
        .collect()
        .filter { bams -> bams }
        .map { bams -> [ [ id: 'consensus_peaks', single_end: false ], bams.toSorted { it.name } ] }

    ch_featurecounts_input = ch_se_batch
        .mix(ch_pe_batch)
        .combine(ch_consensus_saf)

    SUBREAD_FEATURECOUNTS (
        ch_featurecounts_input
    )

    //
    // Merge the per-library-type count matrices into a single consensus matrix
    //
    // Sorted for the same reason: the merge script's column order follows the
    // order of the per-batch matrices it is handed.
    ch_merged_counts = SUBREAD_FEATURECOUNTS.out.counts
        .map { _meta, counts -> counts }
        .collect()
        .map { counts -> [ [ id: 'consensus_peaks' ], counts.toSorted { it.name } ] }

    FEATURECOUNTS_MERGE (
        ch_merged_counts
    )

    //
    // Generate QC plots with DESeq2
    //
    ch_deseq2_qc_pdf           = channel.empty()
    ch_deseq2_qc_rdata         = channel.empty()
    ch_deseq2_qc_rds           = channel.empty()
    ch_deseq2_qc_pca_txt       = channel.empty()
    ch_deseq2_qc_pca_multiqc   = channel.empty()
    ch_deseq2_qc_dists_txt     = channel.empty()
    ch_deseq2_qc_dists_multiqc = channel.empty()
    ch_deseq2_qc_log           = channel.empty()
    ch_deseq2_qc_size_factors  = channel.empty()
    if (!skip_deseq2_qc) {
        DESEQ2_QC (
            FEATURECOUNTS_MERGE.out.counts,
            ch_deseq2_pca_header_multiqc,
            ch_deseq2_clustering_header_multiqc
        )
        ch_deseq2_qc_pdf           = DESEQ2_QC.out.pdf
        ch_deseq2_qc_rdata         = DESEQ2_QC.out.rdata
        ch_deseq2_qc_rds           = DESEQ2_QC.out.rds
        ch_deseq2_qc_pca_txt       = DESEQ2_QC.out.pca_txt
        ch_deseq2_qc_pca_multiqc   = DESEQ2_QC.out.pca_multiqc
        ch_deseq2_qc_dists_txt     = DESEQ2_QC.out.dists_txt
        ch_deseq2_qc_dists_multiqc = DESEQ2_QC.out.dists_multiqc
        ch_deseq2_qc_log           = DESEQ2_QC.out.log
        ch_deseq2_qc_size_factors  = DESEQ2_QC.out.size_factors
    }

    emit:
    consensus_bed           = MACS3_CONSENSUS.out.bed           // channel: [ bed ]
    consensus_saf           = MACS3_CONSENSUS.out.saf           // channel: [ saf ]
    consensus_pdf           = MACS3_CONSENSUS.out.pdf           // channel: [ pdf ]
    consensus_boolean_txt   = MACS3_CONSENSUS.out.boolean_txt   // channel: [ txt ]
    consensus_intersect_txt = MACS3_CONSENSUS.out.intersect_txt // channel: [ txt ]

    homer_annotatepeaks     = ch_homer_annotatepeaks            // channel: [ txt ]

    featurecounts_txt       = FEATURECOUNTS_MERGE.out.counts    // channel: [ val(meta), txt ]
    featurecounts_summary   = SUBREAD_FEATURECOUNTS.out.summary // channel: [ val(meta), txt ] (one per library type)

    deseq2_qc_pdf           = ch_deseq2_qc_pdf                  // channel: [ pdf ]
    deseq2_qc_rdata         = ch_deseq2_qc_rdata                // channel: [ rdata ]
    deseq2_qc_rds           = ch_deseq2_qc_rds                  // channel: [ rds ]
    deseq2_qc_pca_txt       = ch_deseq2_qc_pca_txt              // channel: [ txt ]
    deseq2_qc_pca_multiqc   = ch_deseq2_qc_pca_multiqc          // channel: [ txt ]
    deseq2_qc_dists_txt     = ch_deseq2_qc_dists_txt            // channel: [ txt ]
    deseq2_qc_dists_multiqc = ch_deseq2_qc_dists_multiqc        // channel: [ txt ]
    deseq2_qc_log           = ch_deseq2_qc_log                  // channel: [ txt ]
    deseq2_qc_size_factors  = ch_deseq2_qc_size_factors         // channel: [ txt ]

}
