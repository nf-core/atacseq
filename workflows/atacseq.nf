/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { IGV     } from '../modules/local/igv'
include { MULTIQC } from '../modules/local/multiqc'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_atacseq_pipeline'
include { INPUT_CHECK            } from '../subworkflows/local/input_check'
include { ALIGN_STAR             } from '../subworkflows/local/align_star'
include { BAM_SHIFT_READS as MERGED_LIBRARY_BAM_SHIFT_READS                   } from '../subworkflows/local/bam_shift_reads'
include { BAM_SHIFT_READS as MERGED_REPLICATE_BAM_SHIFT_READS                 } from '../subworkflows/local/bam_shift_reads'
include { BIGWIG_PLOT_DEEPTOOLS as MERGED_LIBRARY_BIGWIG_PLOT_DEEPTOOLS       } from '../subworkflows/local/bigwig_plot_deeptools'
include { BAM_FILTER_BAMTOOLS as MERGED_LIBRARY_FILTER_BAM                    } from '../subworkflows/local/bam_filter_bamtools'
include { BAM_BEDGRAPH_BIGWIG_BEDTOOLS_UCSC as MERGED_LIBRARY_BAM_TO_BIGWIG   } from '../subworkflows/local/bam_bedgraph_bigwig_bedtools_ucsc'
include { BAM_BEDGRAPH_BIGWIG_BEDTOOLS_UCSC as MERGED_REPLICATE_BAM_TO_BIGWIG } from '../subworkflows/local/bam_bedgraph_bigwig_bedtools_ucsc'

include { BAM_PEAKS_CALL_QC_ANNOTATE_MACS2_HOMER as MERGED_LIBRARY_CALL_ANNOTATE_PEAKS   } from '../subworkflows/local/bam_peaks_call_qc_annotate_macs2_homer.nf'
include { BAM_PEAKS_CALL_QC_ANNOTATE_MACS2_HOMER as MERGED_REPLICATE_CALL_ANNOTATE_PEAKS } from '../subworkflows/local/bam_peaks_call_qc_annotate_macs2_homer.nf'
include { BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2 as MERGED_LIBRARY_CONSENSUS_PEAKS   } from '../subworkflows/local/bed_consensus_quantify_qc_bedtools_featurecounts_deseq2.nf'
include { BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2 as MERGED_REPLICATE_CONSENSUS_PEAKS } from '../subworkflows/local/bed_consensus_quantify_qc_bedtools_featurecounts_deseq2.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { PICARD_COLLECTMULTIPLEMETRICS as MERGED_LIBRARY_PICARD_COLLECTMULTIPLEMETRICS } from '../modules/nf-core/picard/collectmultiplemetrics/main'
include { PRESEQ_LCEXTRAP as MERGED_LIBRARY_PRESEQ_LCEXTRAP                             } from '../modules/nf-core/preseq/lcextrap/main'
include { DEEPTOOLS_PLOTFINGERPRINT as MERGED_LIBRARY_DEEPTOOLS_PLOTFINGERPRINT         } from '../modules/nf-core/deeptools/plotfingerprint/main'
include { ATAQV_ATAQV as MERGED_LIBRARY_ATAQV_ATAQV                                     } from '../modules/nf-core/ataqv/ataqv/main'
include { ATAQV_MKARV as MERGED_LIBRARY_ATAQV_MKARV                                     } from '../modules/nf-core/ataqv/mkarv/main'

include { PICARD_MERGESAMFILES as PICARD_MERGESAMFILES_LIBRARY   } from '../modules/nf-core/picard/mergesamfiles/main'
include { PICARD_MERGESAMFILES as PICARD_MERGESAMFILES_REPLICATE } from '../modules/nf-core/picard/mergesamfiles/main'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main'
include { FASTQ_ALIGN_BWA                  } from '../subworkflows/nf-core/fastq_align_bwa/main'
include { FASTQ_ALIGN_BOWTIE2              } from '../subworkflows/nf-core/fastq_align_bowtie2/main'
include { FASTQ_ALIGN_CHROMAP              } from '../subworkflows/nf-core/fastq_align_chromap/main'

include { BAM_MARKDUPLICATES_PICARD as MERGED_LIBRARY_MARKDUPLICATES_PICARD   } from '../subworkflows/nf-core/bam_markduplicates_picard/main'
include { BAM_MARKDUPLICATES_PICARD as MERGED_REPLICATE_MARKDUPLICATES_PICARD } from '../subworkflows/nf-core/bam_markduplicates_picard/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// JSON files required by BAMTools for alignment filtering
ch_bamtools_filter_se_config = file(params.bamtools_filter_se_config)
ch_bamtools_filter_pe_config = file(params.bamtools_filter_pe_config)

// Header files for MultiQC
ch_multiqc_merged_library_peak_count_header        = file("$projectDir/assets/multiqc/merged_library_peak_count_header.txt", checkIfExists: true)
ch_multiqc_merged_library_frip_score_header        = file("$projectDir/assets/multiqc/merged_library_frip_score_header.txt", checkIfExists: true)
ch_multiqc_merged_library_peak_annotation_header   = file("$projectDir/assets/multiqc/merged_library_peak_annotation_header.txt", checkIfExists: true)
ch_multiqc_merged_library_deseq2_pca_header        = file("$projectDir/assets/multiqc/merged_library_deseq2_pca_header.txt", checkIfExists: true)
ch_multiqc_merged_library_deseq2_clustering_header = file("$projectDir/assets/multiqc/merged_library_deseq2_clustering_header.txt", checkIfExists: true)

ch_multiqc_merged_replicate_peak_count_header        = file("$projectDir/assets/multiqc/merged_replicate_peak_count_header.txt", checkIfExists: true)
ch_multiqc_merged_replicate_frip_score_header        = file("$projectDir/assets/multiqc/merged_replicate_frip_score_header.txt", checkIfExists: true)
ch_multiqc_merged_replicate_peak_annotation_header   = file("$projectDir/assets/multiqc/merged_replicate_peak_annotation_header.txt", checkIfExists: true)
ch_multiqc_merged_replicate_deseq2_pca_header        = file("$projectDir/assets/multiqc/merged_replicate_deseq2_pca_header.txt", checkIfExists: true)
ch_multiqc_merged_replicate_deseq2_clustering_header = file("$projectDir/assets/multiqc/merged_replicate_deseq2_clustering_header.txt", checkIfExists: true)

// Check ataqv_mito_reference parameter
ataqv_mito_reference = params.ataqv_mito_reference
if (!params.ataqv_mito_reference && params.mito_name) {
    ataqv_mito_reference = params.mito_name
}

workflow ATACSEQ {

    take:
    ch_samplesheet   // channel: path(sample_sheet.csv)
    ch_versions      // channel: [ path(versions.yml) ]
    ch_fasta         // channel: path(genome.fa)
    ch_fai           // channel: path(genome.fai)
    ch_gtf           // channel: path(genome.gtf)
    ch_gene_bed      // channel: path(gene.beds)
    ch_tss_bed       // channel: path(genome.tss.bed)
    ch_chrom_sizes   // channel: path(chrom.sizes)
    ch_filtered_bed  // channel: path(filtered.bed)
    ch_bwa_index     // channel: path(bwa/index/)
    ch_bowtie2_index // channel: path(bowtie2/index)
    ch_chromap_index // channel: path(chromap.index)
    ch_star_index    // channel: path(star/index/)
    ch_autosomes     // channel: path(autosomes.txt)
    ch_macs_gsize    // channel: integer

    main:
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    ch_input = file(ch_samplesheet)
    INPUT_CHECK (
        ch_input,
        params.seq_center,
        params.with_control
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // Check if reads are all paired-end if 'shift_reads' parameter is set
    //
    if (params.shift_reads) {
        INPUT_CHECK
            .out
            .reads
            .filter { meta, reads -> meta.single_end }
            .collect()
            .map {
                it ->
                    def count = it.size()
                    if (count > 0) {
                        exit 1, 'The parameter --shift_reads can only be applied if all samples are paired-end.'
                    }
            }
    }

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQ_FASTQC_UMITOOLS_TRIMGALORE (
        INPUT_CHECK.out.reads,
        params.skip_fastqc || params.skip_qc,
        false,
        false,
        params.skip_trimming,
        0,
        params.min_trimmed_reads
    )
    ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions)

    //
    // SUBWORKFLOW: Alignment with BWA & BAM QC
    //
    ch_genome_bam        = Channel.empty()
    ch_samtools_stats    = Channel.empty()
    ch_samtools_flagstat = Channel.empty()
    ch_samtools_idxstats = Channel.empty()
    if (params.aligner == 'bwa') {
        FASTQ_ALIGN_BWA (
            FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads,
            ch_bwa_index,
            false,
            ch_fasta
                .map {
                        [ [:], it ]
                }
        )
        ch_genome_bam        = FASTQ_ALIGN_BWA.out.bam
        ch_samtools_stats    = FASTQ_ALIGN_BWA.out.stats
        ch_samtools_flagstat = FASTQ_ALIGN_BWA.out.flagstat
        ch_samtools_idxstats = FASTQ_ALIGN_BWA.out.idxstats
        ch_versions = ch_versions.mix(FASTQ_ALIGN_BWA.out.versions)
    }

    //
    // SUBWORKFLOW: Alignment with BOWTIE2 & BAM QC
    //
    if (params.aligner == 'bowtie2') {
        FASTQ_ALIGN_BOWTIE2 (
            FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads,
            ch_bowtie2_index,
            params.save_unaligned,
            false,
            ch_fasta
                .map {
                    [ [:], it ]
                }
        )
        ch_genome_bam        = FASTQ_ALIGN_BOWTIE2.out.bam
        ch_samtools_stats    = FASTQ_ALIGN_BOWTIE2.out.stats
        ch_samtools_flagstat = FASTQ_ALIGN_BOWTIE2.out.flagstat
        ch_samtools_idxstats = FASTQ_ALIGN_BOWTIE2.out.idxstats
        ch_versions = ch_versions.mix(FASTQ_ALIGN_BOWTIE2.out.versions)
    }

    //
    // SUBWORKFLOW: Alignment with CHROMAP & BAM QC
    //
    if (params.aligner == 'chromap') {
        FASTQ_ALIGN_CHROMAP (
            FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads,
            ch_chromap_index,
            ch_fasta
                .map {
                    [ [:], it ]
                },
            [],
            [],
            [],
            []
        )
        ch_genome_bam        = FASTQ_ALIGN_CHROMAP.out.bam
        ch_genome_bam_index  = FASTQ_ALIGN_CHROMAP.out.bai
        ch_samtools_stats    = FASTQ_ALIGN_CHROMAP.out.stats
        ch_samtools_flagstat = FASTQ_ALIGN_CHROMAP.out.flagstat
        ch_samtools_idxstats = FASTQ_ALIGN_CHROMAP.out.idxstats
        ch_versions = ch_versions.mix(FASTQ_ALIGN_CHROMAP.out.versions)
    }

    //
    // SUBWORKFLOW: Alignment with STAR & BAM QC
    //
    ch_star_multiqc = Channel.empty()
    if (params.aligner == 'star') {
        ALIGN_STAR (
            FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads,
            ch_star_index,
            ch_fasta
                .map {
                    [ [:], it ]
                },
            params.seq_center ?: ''
        )
        ch_genome_bam        = ALIGN_STAR.out.bam
        ch_samtools_stats    = ALIGN_STAR.out.stats
        ch_samtools_flagstat = ALIGN_STAR.out.flagstat
        ch_samtools_idxstats = ALIGN_STAR.out.idxstats
        ch_star_multiqc      = ALIGN_STAR.out.log_final
        ch_versions = ch_versions.mix(ALIGN_STAR.out.versions)
    }

    // Create channels: [ meta, [bam] ]
    ch_genome_bam
        .map {
            meta, bam ->
                def meta_clone = meta.clone()
                meta_clone.remove('read_group')
                meta_clone.id = meta_clone.id - ~/_T\d+$/
                [ meta_clone, bam ]
        }
        .groupTuple(by: [0])
        .map {
            meta, bam ->
                [ meta, bam.flatten() ]
        }
        .set { ch_sort_bam }

    //
    // MODULE: Merge resequenced BAM files
    //
    PICARD_MERGESAMFILES_LIBRARY (
        ch_sort_bam
    )
    ch_versions = ch_versions.mix(PICARD_MERGESAMFILES_LIBRARY.out.versions.first())

    //
    // SUBWORKFLOW: Mark duplicates in BAM files
    //
    MERGED_LIBRARY_MARKDUPLICATES_PICARD (
        PICARD_MERGESAMFILES_LIBRARY.out.bam,
        ch_fasta
            .map {
                [ [:], it ]
            },
        ch_fai
            .map {
                [ [:], it ]
            }
    )
    ch_versions = ch_versions.mix(MERGED_LIBRARY_MARKDUPLICATES_PICARD.out.versions)

    //
    // SUBWORKFLOW: Filter BAM file
    //
    MERGED_LIBRARY_FILTER_BAM (
        MERGED_LIBRARY_MARKDUPLICATES_PICARD.out.bam
            .join(MERGED_LIBRARY_MARKDUPLICATES_PICARD.out.bai, by: [0], remainder: true)
            .join(MERGED_LIBRARY_MARKDUPLICATES_PICARD.out.csi, by: [0], remainder: true)
            .map {
                meta, bam, bai, csi ->
                    if (bai) {
                        [ meta, bam, bai ]
                    } else {
                        [ meta, bam, csi ]
                    }
            },
        ch_filtered_bed.first(),
        ch_fasta
            .map {
                [ [:], it ]
            },
        ch_bamtools_filter_se_config,
        ch_bamtools_filter_pe_config
    )
    ch_versions = ch_versions.mix(MERGED_LIBRARY_FILTER_BAM.out.versions)

    //
    // MODULE: Preseq coverage analysis
    //
    ch_preseq_multiqc = Channel.empty()
    if (!params.skip_preseq) {
        MERGED_LIBRARY_PRESEQ_LCEXTRAP (
            MERGED_LIBRARY_MARKDUPLICATES_PICARD.out.bam
        )
        ch_preseq_multiqc = MERGED_LIBRARY_PRESEQ_LCEXTRAP.out.lc_extrap
        ch_versions = ch_versions.mix(MERGED_LIBRARY_PRESEQ_LCEXTRAP.out.versions.first())
    }

    //
    // MODULE: Picard post alignment QC
    //
    ch_picardcollectmultiplemetrics_multiqc = Channel.empty()
    if (!params.skip_picard_metrics) {
        MERGED_LIBRARY_PICARD_COLLECTMULTIPLEMETRICS (
            MERGED_LIBRARY_FILTER_BAM
                .out
                .bam
                .map {
                    [ it[0], it[1], [] ]
                },
            ch_fasta
                .map {
                    [ [:], it ]
                },
            ch_fai
                .map {
                    [ [:], it ]
                }
        )
        ch_picardcollectmultiplemetrics_multiqc = MERGED_LIBRARY_PICARD_COLLECTMULTIPLEMETRICS.out.metrics
        ch_versions = ch_versions.mix(MERGED_LIBRARY_PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())
    }

    //
    // SUBWORKFLOW: Shift paired-end reads
    //
    ch_merged_library_filter_bam      = MERGED_LIBRARY_FILTER_BAM.out.bam
    ch_merged_library_filter_bai      = MERGED_LIBRARY_FILTER_BAM.out.bai
    ch_merged_library_filter_flagstat = MERGED_LIBRARY_FILTER_BAM.out.flagstat
    ch_merged_library_filter_csi      = MERGED_LIBRARY_FILTER_BAM.out.csi

    if (params.shift_reads && params.aligner != 'chromap' ) {
        MERGED_LIBRARY_BAM_SHIFT_READS (
            ch_merged_library_filter_bam.join(ch_merged_library_filter_bai, by: [0]),
            ch_fasta
            .map {
                [ [:], it ]
            }
        )
        ch_versions = ch_versions.mix(MERGED_LIBRARY_BAM_SHIFT_READS.out.versions)

        ch_merged_library_filter_bam      = MERGED_LIBRARY_BAM_SHIFT_READS.out.bam
        ch_merged_library_filter_bai      = MERGED_LIBRARY_BAM_SHIFT_READS.out.bai
        ch_merged_library_filter_flagstat = MERGED_LIBRARY_BAM_SHIFT_READS.out.flagstat
        ch_merged_library_filter_csi      = MERGED_LIBRARY_BAM_SHIFT_READS.out.csi

    }

    //
    // SUBWORKFLOW: Normalised bigWig coverage tracks
    //
    MERGED_LIBRARY_BAM_TO_BIGWIG (
        // MERGED_LIBRARY_FILTER_BAM.out.bam.join(MERGED_LIBRARY_FILTER_BAM.out.flagstat, by: [0]),
        // ch_chrom_sizes
        ch_merged_library_filter_bam.join(ch_merged_library_filter_flagstat, by: [0]),
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(MERGED_LIBRARY_BAM_TO_BIGWIG.out.versions)

    //
    // SUBWORKFLOW: Plot coverage across annotation with deepTools
    //
    ch_deeptoolsplotprofile_multiqc = Channel.empty()
    if (!params.skip_plot_profile) {
        MERGED_LIBRARY_BIGWIG_PLOT_DEEPTOOLS (
            MERGED_LIBRARY_BAM_TO_BIGWIG.out.bigwig,
            ch_gene_bed,
            ch_tss_bed
        )
        ch_deeptoolsplotprofile_multiqc = MERGED_LIBRARY_BIGWIG_PLOT_DEEPTOOLS.out.plotprofile_table
        ch_versions = ch_versions.mix(MERGED_LIBRARY_BIGWIG_PLOT_DEEPTOOLS.out.versions)
    }

    // Create channels: [ meta, [bam], [bai] ] or [ meta, [ bam, control_bam ] [ bai, control_bai ] ]
    ch_merged_library_filter_bam
        .join(ch_merged_library_filter_bai, by: [0], remainder: true)
        .join(ch_merged_library_filter_csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }
        .set { ch_bam_bai }

    if (params.with_control) {
        ch_bam_bai
            .map {
                meta, bam, bai ->
                    meta.control ? null : [ meta.id, [ bam ] , [ bai ] ]
            }
            .set { ch_control_bam_bai }

        ch_bam_bai
            .map {
                meta, bam, bai ->
                    meta.control ? [ meta.control, meta, [ bam ], [ bai ] ] : null
            }
            .combine(ch_control_bam_bai, by: 0)
            .map { it -> [ it[1] , it[2] + it[4], it[3] + it[5] ] }
            .set { ch_bam_bai }
    }

    //
    // MODULE: deepTools plotFingerprint QC
    //
    ch_deeptoolsplotfingerprint_multiqc = Channel.empty()
    if (!params.skip_plot_fingerprint) {
        MERGED_LIBRARY_DEEPTOOLS_PLOTFINGERPRINT (
            ch_bam_bai
        )
        ch_deeptoolsplotfingerprint_multiqc = MERGED_LIBRARY_DEEPTOOLS_PLOTFINGERPRINT.out.matrix
        ch_versions = ch_versions.mix(MERGED_LIBRARY_DEEPTOOLS_PLOTFINGERPRINT.out.versions.first())
    }

    // Create channel: [ val(meta), bam, control_bam ]
    if (params.with_control) {
        ch_bam_bai
            .map {
                meta, bams, bais ->
                    [ meta , bams[0], bams[1] ]
            }
            .set { ch_bam_library }
    } else {
        ch_bam_bai
            .map {
                meta, bam, bai ->
                    [ meta , bam, [] ]
            }
            .set { ch_bam_library }
    }

    //
    // SUBWORKFLOW: Call peaks with MACS2, annotate with HOMER and perform downstream QC
    //
    MERGED_LIBRARY_CALL_ANNOTATE_PEAKS (
        ch_bam_library,
        ch_fasta,
        ch_gtf,
        ch_macs_gsize,
        ".mLb.clN_peaks.annotatePeaks.txt",
        ch_multiqc_merged_library_peak_count_header,
        ch_multiqc_merged_library_frip_score_header,
        ch_multiqc_merged_library_peak_annotation_header,
        params.narrow_peak,
        params.skip_peak_annotation,
        params.skip_peak_qc
    )
    ch_versions = ch_versions.mix(MERGED_LIBRARY_CALL_ANNOTATE_PEAKS.out.versions)

    //
    // SUBWORKFLOW: Consensus peaks analysis
    //
    ch_macs2_consensus_library_bed       = Channel.empty()
    ch_featurecounts_library_multiqc     = Channel.empty()
    ch_deseq2_pca_library_multiqc        = Channel.empty()
    ch_deseq2_clustering_library_multiqc = Channel.empty()
    if (!params.skip_consensus_peaks) {
        MERGED_LIBRARY_CONSENSUS_PEAKS (
            MERGED_LIBRARY_CALL_ANNOTATE_PEAKS.out.peaks,
            ch_bam_library,
            ch_fasta,
            ch_gtf,
            ch_multiqc_merged_library_deseq2_pca_header,
            ch_multiqc_merged_library_deseq2_clustering_header,
            params.narrow_peak,
            params.skip_peak_annotation,
            params.skip_deseq2_qc
        )
        ch_macs2_consensus_library_bed       = MERGED_LIBRARY_CONSENSUS_PEAKS.out.consensus_bed
        ch_featurecounts_library_multiqc     = MERGED_LIBRARY_CONSENSUS_PEAKS.out.featurecounts_summary
        ch_deseq2_pca_library_multiqc        = MERGED_LIBRARY_CONSENSUS_PEAKS.out.deseq2_qc_pca_multiqc
        ch_deseq2_clustering_library_multiqc = MERGED_LIBRARY_CONSENSUS_PEAKS.out.deseq2_qc_dists_multiqc
        ch_versions = ch_versions.mix(MERGED_LIBRARY_CONSENSUS_PEAKS.out.versions)
    }

    // Create channels: [ meta, bam, bai, peak_file ]
    MERGED_LIBRARY_MARKDUPLICATES_PICARD
        .out
        .bam
        .join(MERGED_LIBRARY_MARKDUPLICATES_PICARD.out.bai, by: [0], remainder: true)
        .join(MERGED_LIBRARY_MARKDUPLICATES_PICARD.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }
        .join(MERGED_LIBRARY_CALL_ANNOTATE_PEAKS.out.peaks, by: [0])
        .set { ch_bam_peaks }

    //
    // MODULE: ATAQV QC
    //
    if (!params.skip_ataqv) {
        MERGED_LIBRARY_ATAQV_ATAQV (
            ch_bam_peaks,
            'NA',
            ataqv_mito_reference ?: '',
            ch_tss_bed,
            [],
            ch_autosomes
        )
        ch_versions = ch_versions.mix(MERGED_LIBRARY_ATAQV_ATAQV.out.versions.first())

        MERGED_LIBRARY_ATAQV_MKARV (
            MERGED_LIBRARY_ATAQV_ATAQV.out.json.collect{it[1]}
        )
        ch_versions = ch_versions.mix(MERGED_LIBRARY_ATAQV_MKARV.out.versions)
    }

    //
    // Merged replicate analysis
    //
    ch_markduplicates_replicate_stats                   = Channel.empty()
    ch_markduplicates_replicate_flagstat                = Channel.empty()
    ch_markduplicates_replicate_idxstats                = Channel.empty()
    ch_markduplicates_replicate_metrics                 = Channel.empty()
    ch_ucsc_bedgraphtobigwig_replicate_bigwig           = Channel.empty()
    ch_macs2_replicate_peaks                            = Channel.empty()
    ch_macs2_frip_replicate_multiqc                     = Channel.empty()
    ch_macs2_peak_count_replicate_multiqc               = Channel.empty()
    ch_macs2_plot_homer_annotatepeaks_replicate_multiqc = Channel.empty()
    ch_macs2_consensus_replicate_bed                    = Channel.empty()
    ch_featurecounts_replicate_multiqc                  = Channel.empty()
    ch_deseq2_pca_replicate_multiqc                     = Channel.empty()
    ch_deseq2_clustering_replicate_multiqc              = Channel.empty()
    if (!params.skip_merge_replicates) {

        // Check if we have multiple replicates
        MERGED_LIBRARY_FILTER_BAM
            .out
            .bam
            .map {
                meta, bam ->
                    def meta_clone = meta.clone()
                    meta_clone.id = meta_clone.id - ~/_REP\d+$/
                    meta_clone.control = meta_clone.control ? meta_clone.control - ~/_REP\d+$/ : ""
                    [ meta_clone.id, meta_clone, bam ]
            }
            .groupTuple()
            .map {
                id, metas, bams ->
                    if (bams.size() > 1) {
                        return [ metas[0], bams ]
                    }
            }
            .set { ch_merged_library_replicate_bam }

        //
        // MODULE: Merge replicate BAM files
        //
        PICARD_MERGESAMFILES_REPLICATE (
            ch_merged_library_replicate_bam
        )
        ch_versions = ch_versions.mix(PICARD_MERGESAMFILES_REPLICATE.out.versions.first())

        //
        // SUBWORKFLOW: Mark duplicates & filter BAM files after merging
        //
        MERGED_REPLICATE_MARKDUPLICATES_PICARD (
            PICARD_MERGESAMFILES_REPLICATE.out.bam,
            ch_fasta
                .map {
                    [ [:], it ]
                },
            ch_fai
                .map {
                    [ [:], it ]
                }
        )
        ch_markduplicates_replicate_stats    = MERGED_REPLICATE_MARKDUPLICATES_PICARD.out.stats
        ch_markduplicates_replicate_flagstat = MERGED_REPLICATE_MARKDUPLICATES_PICARD.out.flagstat
        ch_markduplicates_replicate_idxstats = MERGED_REPLICATE_MARKDUPLICATES_PICARD.out.idxstats
        ch_markduplicates_replicate_metrics  = MERGED_REPLICATE_MARKDUPLICATES_PICARD.out.metrics
        ch_versions = ch_versions.mix(MERGED_REPLICATE_MARKDUPLICATES_PICARD.out.versions)

        //
        // SUBWORKFLOW: Shift paired-end reads
        // Shift again, as ch_merged_library_replicate_bam is generated out of unshifted reads
        //
        ch_merged_replicate_markduplicate_bam = MERGED_REPLICATE_MARKDUPLICATES_PICARD.out.bam
        ch_merged_replicate_markduplicate_bai = MERGED_REPLICATE_MARKDUPLICATES_PICARD.out.bai
        ch_merged_replicate_markduplicate_flagstat = MERGED_REPLICATE_MARKDUPLICATES_PICARD.out.flagstat
        ch_merged_replicate_markduplicate_csi = MERGED_REPLICATE_MARKDUPLICATES_PICARD.out.csi

        if (params.shift_reads && params.aligner != 'chromap' ) {
            MERGED_REPLICATE_BAM_SHIFT_READS (
                ch_merged_replicate_markduplicate_bam.join(ch_merged_replicate_markduplicate_bai, by: [0]),
                ch_fasta
                .map {
                    [ [:], it ]
                }
            )
            ch_versions = ch_versions.mix(MERGED_REPLICATE_BAM_SHIFT_READS.out.versions)

            ch_merged_replicate_markduplicate_bam      = MERGED_REPLICATE_BAM_SHIFT_READS.out.bam
            ch_merged_replicate_markduplicate_bai      = MERGED_REPLICATE_BAM_SHIFT_READS.out.bai
            ch_merged_replicate_markduplicate_flagstat = MERGED_REPLICATE_BAM_SHIFT_READS.out.flagstat
            ch_merged_replicate_markduplicate_csi      = MERGED_REPLICATE_BAM_SHIFT_READS.out.csi
        }
        if (!params.skip_merged_replicate_bigwig) {
            //
            // SUBWORKFLOW: Normalised bigWig coverage tracks
            //
            MERGED_REPLICATE_BAM_TO_BIGWIG (
                ch_merged_replicate_markduplicate_bam.join(ch_merged_replicate_markduplicate_flagstat, by: [0]),
                ch_chrom_sizes
            )
            ch_ucsc_bedgraphtobigwig_replicate_bigwig = MERGED_REPLICATE_BAM_TO_BIGWIG.out.bigwig
            ch_versions = ch_versions.mix(MERGED_REPLICATE_BAM_TO_BIGWIG.out.versions)
        }
        // Create channels: [ meta, bam, ([] for control_bam) ]
        if (params.with_control) {
            ch_merged_replicate_markduplicate_bam
                .map {
                    meta, bam ->
                        meta.control ? null : [ meta.id, bam ]
                }
                .set { ch_bam_merged_control }

            ch_merged_replicate_markduplicate_bam
                .map {
                    meta, bam ->
                        meta.control ? [ meta.control, meta, bam ] : null
                }
                .combine( ch_bam_merged_control, by: 0)
                .map { it -> [ it[1] , it[2], it[3] ] }
                .set { ch_bam_replicate }
        } else {
            ch_merged_replicate_markduplicate_bam
                .map {
                    meta, bam ->
                        [ meta , bam, [] ]
                }
                .set { ch_bam_replicate }
        }

        //
        // SUBWORKFLOW: Call peaks with MACS2, annotate with HOMER and perform downstream QC
        //
        MERGED_REPLICATE_CALL_ANNOTATE_PEAKS (
            ch_bam_replicate,
            ch_fasta,
            ch_gtf,
            ch_macs_gsize,
            ".mRp.clN_peaks.annotatePeaks.txt",
            ch_multiqc_merged_replicate_peak_count_header,
            ch_multiqc_merged_replicate_frip_score_header,
            ch_multiqc_merged_replicate_peak_annotation_header,
            params.narrow_peak,
            params.skip_peak_annotation,
            params.skip_peak_qc
        )
        ch_macs2_replicate_peaks                            = MERGED_REPLICATE_CALL_ANNOTATE_PEAKS.out.peaks
        ch_macs2_frip_replicate_multiqc                     = MERGED_REPLICATE_CALL_ANNOTATE_PEAKS.out.frip_multiqc
        ch_macs2_peak_count_replicate_multiqc               = MERGED_REPLICATE_CALL_ANNOTATE_PEAKS.out.peak_count_multiqc
        ch_macs2_plot_homer_annotatepeaks_replicate_multiqc = MERGED_REPLICATE_CALL_ANNOTATE_PEAKS.out.plot_homer_annotatepeaks_tsv
        ch_versions = ch_versions.mix(MERGED_REPLICATE_CALL_ANNOTATE_PEAKS.out.versions)

        //
        // SUBWORKFLOW: Consensus peaks analysis
        //
        if (!params.skip_consensus_peaks) {
            MERGED_REPLICATE_CONSENSUS_PEAKS (
                MERGED_REPLICATE_CALL_ANNOTATE_PEAKS.out.peaks,
                ch_merged_library_replicate_bam,
                ch_fasta,
                ch_gtf,
                ch_multiqc_merged_replicate_deseq2_pca_header,
                ch_multiqc_merged_replicate_deseq2_clustering_header,
                params.narrow_peak,
                params.skip_peak_annotation,
                params.skip_deseq2_qc
            )
            ch_macs2_consensus_replicate_bed       = MERGED_REPLICATE_CONSENSUS_PEAKS.out.consensus_bed
            ch_featurecounts_replicate_multiqc     = MERGED_REPLICATE_CONSENSUS_PEAKS.out.featurecounts_summary
            ch_deseq2_pca_replicate_multiqc        = MERGED_REPLICATE_CONSENSUS_PEAKS.out.deseq2_qc_pca_multiqc
            ch_deseq2_clustering_replicate_multiqc = MERGED_REPLICATE_CONSENSUS_PEAKS.out.deseq2_qc_dists_multiqc
            ch_versions = ch_versions.mix(MERGED_REPLICATE_CONSENSUS_PEAKS.out.versions)
        }
    }

    //
    // MODULE: Create IGV session
    //
    if (!params.skip_igv) {
        IGV (
            ch_fasta,
            ch_fai,
            MERGED_LIBRARY_BAM_TO_BIGWIG.out.bigwig.collect{it[1]}.ifEmpty([]),
            MERGED_LIBRARY_CALL_ANNOTATE_PEAKS.out.peaks.collect{it[1]}.ifEmpty([]),
            ch_macs2_consensus_library_bed.collect{it[1]}.ifEmpty([]),
            ch_ucsc_bedgraphtobigwig_replicate_bigwig.collect{it[1]}.ifEmpty([]),
            ch_macs2_replicate_peaks.collect{it[1]}.ifEmpty([]),
            ch_macs2_consensus_replicate_bed.collect{it[1]}.ifEmpty([]),
            "${params.aligner}/merged_library/bigwig",
            { ["${params.aligner}/merged_library/macs2",
                params.narrow_peak? '/narrow_peak' : '/broad_peak'
                ].join('') },
            { ["${params.aligner}/merged_library/macs2",
                params.narrow_peak? '/narrow_peak' : '/broad_peak',
                "/consensus"
                ].join('') },
            "${params.aligner}/merged_replicate/bigwig",
            { ["${params.aligner}/merged_replicate/macs2",
                params.narrow_peak? '/narrow_peak' : '/broad_peak'
                ].join('') },
            { ["${params.aligner}/merged_replicate/macs2",
                params.narrow_peak? '/narrow_peak' : '/broad_peak',
                "/consensus"
                ].join('') },
        )
        ch_versions = ch_versions.mix(IGV.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_atacseq_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
        ch_multiqc_logo                       = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo)   : Channel.empty()
        summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
        ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),

            FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
            FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]),
            FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]),

            ch_samtools_stats.collect{it[1]}.ifEmpty([]),
            ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
            ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),

            MERGED_LIBRARY_MARKDUPLICATES_PICARD.out.stats.collect{it[1]}.ifEmpty([]),
            MERGED_LIBRARY_MARKDUPLICATES_PICARD.out.flagstat.collect{it[1]}.ifEmpty([]),
            MERGED_LIBRARY_MARKDUPLICATES_PICARD.out.idxstats.collect{it[1]}.ifEmpty([]),
            MERGED_LIBRARY_MARKDUPLICATES_PICARD.out.metrics.collect{it[1]}.ifEmpty([]),

            MERGED_LIBRARY_FILTER_BAM.out.stats.collect{it[1]}.ifEmpty([]),
            MERGED_LIBRARY_FILTER_BAM.out.flagstat.collect{it[1]}.ifEmpty([]),
            MERGED_LIBRARY_FILTER_BAM.out.idxstats.collect{it[1]}.ifEmpty([]),
            ch_picardcollectmultiplemetrics_multiqc.collect{it[1]}.ifEmpty([]),

            ch_preseq_multiqc.collect{it[1]}.ifEmpty([]),

            ch_deeptoolsplotprofile_multiqc.collect{it[1]}.ifEmpty([]),
            ch_deeptoolsplotfingerprint_multiqc.collect{it[1]}.ifEmpty([]),

            MERGED_LIBRARY_CALL_ANNOTATE_PEAKS.out.frip_multiqc.collect{it[1]}.ifEmpty([]),
            MERGED_LIBRARY_CALL_ANNOTATE_PEAKS.out.peak_count_multiqc.collect{it[1]}.ifEmpty([]),
            MERGED_LIBRARY_CALL_ANNOTATE_PEAKS.out.plot_homer_annotatepeaks_tsv.collect().ifEmpty([]),
            ch_featurecounts_library_multiqc.collect{it[1]}.ifEmpty([]),

            ch_markduplicates_replicate_stats.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_replicate_flagstat.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_replicate_idxstats.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_replicate_metrics.collect{it[1]}.ifEmpty([]),

            ch_macs2_frip_replicate_multiqc.collect{it[1]}.ifEmpty([]),
            ch_macs2_peak_count_replicate_multiqc.collect{it[1]}.ifEmpty([]),
            ch_macs2_plot_homer_annotatepeaks_replicate_multiqc.collect().ifEmpty([]),
            ch_featurecounts_replicate_multiqc.collect{it[1]}.ifEmpty([]),

            ch_deseq2_pca_library_multiqc.collect().ifEmpty([]),
            ch_deseq2_clustering_library_multiqc.collect().ifEmpty([]),
            ch_deseq2_pca_replicate_multiqc.collect().ifEmpty([]),
            ch_deseq2_clustering_replicate_multiqc.collect().ifEmpty([])
        )
        ch_multiqc_report = MULTIQC.out.report
    }

    emit:
    multiqc_report = ch_multiqc_report  // channel: /path/to/multiqc_report.html
    versions       = ch_versions       // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
