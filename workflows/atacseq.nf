/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def valid_params = [
    aligners : [ 'bwa', 'bowtie2', 'chromap', 'star' ]
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowAtacseq.initialise(params, log, valid_params)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input, params.multiqc_config,
    params.fasta,
    params.gtf, params.gff, params.gene_bed,
    params.bwa_index, params.bowtie2_index, params.chromap_index, params.star_index,
    params.blacklist,
    params.bamtools_filter_pe_config, params.bamtools_filter_se_config
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo          = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)


// JSON files required by BAMTools for alignment filtering
ch_bamtools_filter_se_config = file(params.bamtools_filter_se_config, checkIfExists: true)
ch_bamtools_filter_pe_config = file(params.bamtools_filter_pe_config, checkIfExists: true)

// Header files for MultiQC
ch_merged_library_peak_count_header        = file("$projectDir/assets/multiqc/merged_library_peak_count_header.txt", checkIfExists: true)
ch_merged_library_frip_score_header        = file("$projectDir/assets/multiqc/merged_library_frip_score_header.txt", checkIfExists: true)
ch_merged_library_peak_annotation_header   = file("$projectDir/assets/multiqc/merged_library_peak_annotation_header.txt", checkIfExists: true)
ch_merged_library_deseq2_pca_header        = file("$projectDir/assets/multiqc/merged_library_deseq2_pca_header.txt", checkIfExists: true)
ch_merged_library_deseq2_clustering_header = file("$projectDir/assets/multiqc/merged_library_deseq2_clustering_header.txt", checkIfExists: true)

ch_merged_replicate_peak_count_header        = file("$projectDir/assets/multiqc/merged_replicate_peak_count_header.txt", checkIfExists: true)
ch_merged_replicate_frip_score_header        = file("$projectDir/assets/multiqc/merged_replicate_frip_score_header.txt", checkIfExists: true)
ch_merged_replicate_peak_annotation_header   = file("$projectDir/assets/multiqc/merged_replicate_peak_annotation_header.txt", checkIfExists: true)
ch_merged_replicate_deseq2_pca_header        = file("$projectDir/assets/multiqc/merged_replicate_deseq2_pca_header.txt", checkIfExists: true)
ch_merged_replicate_deseq2_clustering_header = file("$projectDir/assets/multiqc/merged_replicate_deseq2_clustering_header.txt", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { IGV      } from '../modules/local/igv'
include { MULTIQC  } from '../modules/local/multiqc'

include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_LIBRARY               } from '../modules/local/bedtools_genomecov'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_REPLICATE             } from '../modules/local/bedtools_genomecov'
include { FRIP_SCORE as FRIP_SCORE_LIBRARY                               } from '../modules/local/frip_score'
include { FRIP_SCORE as FRIP_SCORE_REPLICATE                             } from '../modules/local/frip_score'
include { PLOT_MACS2_QC as PLOT_MACS2_QC_LIBRARY                         } from '../modules/local/plot_macs2_qc'
include { PLOT_MACS2_QC as PLOT_MACS2_QC_REPLICATE                       } from '../modules/local/plot_macs2_qc'
include { PLOT_HOMER_ANNOTATEPEAKS as PLOT_HOMER_ANNOTATEPEAKS_LIBRARY   } from '../modules/local/plot_homer_annotatepeaks'
include { PLOT_HOMER_ANNOTATEPEAKS as PLOT_HOMER_ANNOTATEPEAKS_REPLICATE } from '../modules/local/plot_homer_annotatepeaks'
include { MACS2_CONSENSUS as MACS2_CONSENSUS_LIBRARY                     } from '../modules/local/macs2_consensus'
include { MACS2_CONSENSUS as MACS2_CONSENSUS_REPLICATE                   } from '../modules/local/macs2_consensus'
include { DESEQ2_QC as DESEQ2_QC_LIBRARY                                 } from '../modules/local/deseq2_qc'
include { DESEQ2_QC as DESEQ2_QC_REPLICATE                               } from '../modules/local/deseq2_qc'
include { MULTIQC_CUSTOM_PEAKS as MULTIQC_CUSTOM_PEAKS_LIBRARY           } from '../modules/local/multiqc_custom_peaks'
include { MULTIQC_CUSTOM_PEAKS as MULTIQC_CUSTOM_PEAKS_REPLICATE         } from '../modules/local/multiqc_custom_peaks'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK         } from '../subworkflows/local/input_check'
include { PREPARE_GENOME      } from '../subworkflows/local/prepare_genome'
include { FILTER_BAM_BAMTOOLS } from '../subworkflows/local/filter_bam_bamtools'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { PICARD_COLLECTMULTIPLEMETRICS } from '../modules/nf-core/picard/collectmultiplemetrics/main'
include { PRESEQ_LCEXTRAP               } from '../modules/nf-core/preseq/lcextrap/main'
include { DEEPTOOLS_PLOTPROFILE         } from '../modules/nf-core/deeptools/plotprofile/main'
include { DEEPTOOLS_PLOTHEATMAP         } from '../modules/nf-core/deeptools/plotheatmap/main'
include { DEEPTOOLS_PLOTFINGERPRINT     } from '../modules/nf-core/deeptools/plotfingerprint/main'
include { KHMER_UNIQUEKMERS             } from '../modules/nf-core/khmer/uniquekmers/main'
include { ATAQV_ATAQV                   } from '../modules/nf-core/ataqv/ataqv/main'
include { ATAQV_MKARV                   } from '../modules/nf-core/ataqv/mkarv/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/nf-core/custom/dumpsoftwareversions/main'


include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_SCALE_REGIONS   } from '../modules/nf-core/deeptools/computematrix/main'
include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_REFERENCE_POINT } from '../modules/nf-core/deeptools/computematrix/main'
include { PICARD_MERGESAMFILES as PICARD_MERGESAMFILES_LIBRARY               } from '../modules/nf-core/picard/mergesamfiles/main'
include { PICARD_MERGESAMFILES as PICARD_MERGESAMFILES_REPLICATE             } from '../modules/nf-core/picard/mergesamfiles/main'
include { UCSC_BEDGRAPHTOBIGWIG as UCSC_BEDGRAPHTOBIGWIG_LIBRARY             } from '../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { UCSC_BEDGRAPHTOBIGWIG as UCSC_BEDGRAPHTOBIGWIG_REPLICATE           } from '../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { MACS2_CALLPEAK as MACS2_CALLPEAK_LIBRARY                           } from '../modules/nf-core/macs2/callpeak/main'
include { MACS2_CALLPEAK as MACS2_CALLPEAK_REPLICATE                         } from '../modules/nf-core/macs2/callpeak/main'
include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_MACS2_LIBRARY           } from '../modules/nf-core/homer/annotatepeaks/main'
include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_MACS2_REPLICATE         } from '../modules/nf-core/homer/annotatepeaks/main'
include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_CONSENSUS_LIBRARY       } from '../modules/nf-core/homer/annotatepeaks/main'
include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_CONSENSUS_REPLICATE     } from '../modules/nf-core/homer/annotatepeaks/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_LIBRARY             } from '../modules/nf-core/subread/featurecounts/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_REPLICATE           } from '../modules/nf-core/subread/featurecounts/main'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//

include { FASTQC_TRIMGALORE } from '../subworkflows/nf-core/fastqc_trimgalore'
include { ALIGN_BWA_MEM     } from '../subworkflows/nf-core/align_bwa_mem'
include { ALIGN_BOWTIE2     } from '../subworkflows/nf-core/align_bowtie2'
include { ALIGN_CHROMAP     } from '../subworkflows/nf-core/align_chromap'
include { ALIGN_STAR        } from '../subworkflows/nf-core/align_star'

include { MARK_DUPLICATES_PICARD as MARK_DUPLICATES_PICARD_LIBRARY   } from '../subworkflows/nf-core/mark_duplicates_picard'
include { MARK_DUPLICATES_PICARD as MARK_DUPLICATES_PICARD_REPLICATE } from '../subworkflows/nf-core/mark_duplicates_picard'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow ATACSEQ {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    PREPARE_GENOME (
        params.aligner
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input,
        params.seq_center
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQC_TRIMGALORE (
        INPUT_CHECK.out.reads,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)

    //
    // SUBWORKFLOW: Alignment with BWA & BAM QC
    //
    ch_genome_bam        = Channel.empty()
    ch_genome_bam_index  = Channel.empty()
    ch_samtools_stats    = Channel.empty()
    ch_samtools_flagstat = Channel.empty()
    ch_samtools_idxstats = Channel.empty()
    if (params.aligner == 'bwa') {
        ALIGN_BWA_MEM (
            FASTQC_TRIMGALORE.out.reads,
            PREPARE_GENOME.out.bwa_index
        )
        ch_genome_bam        = ALIGN_BWA_MEM.out.bam
        ch_genome_bam_index  = ALIGN_BWA_MEM.out.bai
        ch_samtools_stats    = ALIGN_BWA_MEM.out.stats
        ch_samtools_flagstat = ALIGN_BWA_MEM.out.flagstat
        ch_samtools_idxstats = ALIGN_BWA_MEM.out.idxstats
        ch_versions = ch_versions.mix(ALIGN_BWA_MEM.out.versions.first())
    }

    //
    // SUBWORKFLOW: Alignment with BOWTIE2 & BAM QC
    //
    if (params.aligner == 'bowtie2') {
        ALIGN_BOWTIE2 (
            FASTQC_TRIMGALORE.out.reads,
            PREPARE_GENOME.out.bowtie2_index,
            params.save_unaligned
        )
        ch_genome_bam        = ALIGN_BOWTIE2.out.bam
        ch_genome_bam_index  = ALIGN_BOWTIE2.out.bai
        ch_samtools_stats    = ALIGN_BOWTIE2.out.stats
        ch_samtools_flagstat = ALIGN_BOWTIE2.out.flagstat
        ch_samtools_idxstats = ALIGN_BOWTIE2.out.idxstats
        ch_versions = ch_versions.mix(ALIGN_BOWTIE2.out.versions.first())
    }

    //
    // SUBWORKFLOW: Alignment with CHROMAP & BAM QC
    //
    if (params.aligner == 'chromap') {
        ALIGN_CHROMAP (
            FASTQC_TRIMGALORE.out.reads,
            PREPARE_GENOME.out.chromap_index,
            PREPARE_GENOME.out.fasta
        )

        // Filter out paired-end reads until the issue below is fixed
        // https://github.com/nf-core/chipseq/issues/291
        // ch_genome_bam = ALIGN_CHROMAP.out.bam
        ALIGN_CHROMAP
            .out
            .bam
            .branch {
                meta, bam ->
                    single_end: meta.single_end
                        return [ meta, bam ]
                    paired_end: !meta.single_end
                        return [ meta, bam ]
            }
            .set { ch_genome_bam_chromap }

        ch_genome_bam_chromap
            .paired_end
            .collect()
            .map {
                it ->
                    def count = it.size()
                    if (count > 0) {
                        log.warn "=============================================================================\n" +
                        "  Paired-end files produced by chromap cannot be used by some downstream tools due to the issue below:\n" +
                        "  https://github.com/nf-core/chipseq/issues/291\n" +
                        "  They will be excluded from the analysis. Consider using a different aligner\n" +
                        "==================================================================================="
                    }
            }

        ch_genome_bam        = ch_genome_bam_chromap.single_end
        ch_genome_bam_index  = ALIGN_CHROMAP.out.bai
        ch_samtools_stats    = ALIGN_CHROMAP.out.stats
        ch_samtools_flagstat = ALIGN_CHROMAP.out.flagstat
        ch_samtools_idxstats = ALIGN_CHROMAP.out.idxstats
        ch_versions = ch_versions.mix(ALIGN_CHROMAP.out.versions.first())
    }

    //
    // SUBWORKFLOW: Alignment with STAR & BAM QC
    //
    if (params.aligner == 'star') {
        ALIGN_STAR (
            FASTQC_TRIMGALORE.out.reads,
            PREPARE_GENOME.out.star_index
        )
        ch_genome_bam        = ALIGN_STAR.out.bam
        ch_genome_bam_index  = ALIGN_STAR.out.bai
        ch_transcriptome_bam = ALIGN_STAR.out.bam_transcript
        ch_samtools_stats    = ALIGN_STAR.out.stats
        ch_samtools_flagstat = ALIGN_STAR.out.flagstat
        ch_samtools_idxstats = ALIGN_STAR.out.idxstats
        ch_star_multiqc      = ALIGN_STAR.out.log_final
        ch_versions = ch_versions.mix(ALIGN_STAR.out.versions)
    }

    //
    // MERGE LIBRARY BAM
    // MODULE: Merge resequenced BAM files
    //
    ch_genome_bam
        .map {
            meta, bam ->
                def fmeta = meta.findAll { it.key != 'read_group' }
                fmeta.id = fmeta.id.split('_')[0..-2].join('_')
                [ fmeta, bam ] }
        .groupTuple(by: [0])
        .map {
            it ->
                [ it[0], it[1].flatten() ]
        }
        .set { ch_sort_bam }

    PICARD_MERGESAMFILES_LIBRARY (
        ch_sort_bam
    )
    ch_versions = ch_versions.mix(PICARD_MERGESAMFILES_LIBRARY.out.versions.first().ifEmpty(null))

    //
    // SUBWORKFLOW: Mark duplicates & filter BAM files after merging
    //
    MARK_DUPLICATES_PICARD_LIBRARY (
        PICARD_MERGESAMFILES_LIBRARY.out.bam
    )
    ch_versions = ch_versions.mix(MARK_DUPLICATES_PICARD_LIBRARY.out.versions)

    //
    // SUBWORKFLOW: Fix getting name sorted BAM here for PE/SE
    //
    FILTER_BAM_BAMTOOLS (
        MARK_DUPLICATES_PICARD_LIBRARY.out.bam.join(MARK_DUPLICATES_PICARD_LIBRARY.out.bai, by: [0]),
        PREPARE_GENOME.out.filtered_bed.first(),

        ch_bamtools_filter_se_config,
        ch_bamtools_filter_pe_config
    )
    ch_versions = ch_versions.mix(FILTER_BAM_BAMTOOLS.out.versions.first().ifEmpty(null))

    //
    // MODULE: Preseq coverage analysis
    //
    ch_preseq_multiqc = Channel.empty()
    if (!params.skip_preseq) {
        PRESEQ_LCEXTRAP (
            MARK_DUPLICATES_PICARD_LIBRARY.out.bam
        )
        ch_preseq_multiqc = PRESEQ_LCEXTRAP.out.lc_extrap
        ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())
    }

    //
    // MODULE: Picard post alignment QC
    //
    ch_picardcollectmultiplemetrics_multiqc = Channel.empty()
    if (!params.skip_picard_metrics) {
        PICARD_COLLECTMULTIPLEMETRICS (
            FILTER_BAM_BAMTOOLS.out.bam,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.fai,
        )
        ch_picardcollectmultiplemetrics_multiqc = PICARD_COLLECTMULTIPLEMETRICS.out.metrics
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())
    }

    //
    // MODULE: BedGraph coverage tracks
    //
    BEDTOOLS_GENOMECOV_LIBRARY (
        FILTER_BAM_BAMTOOLS.out.bam.join(FILTER_BAM_BAMTOOLS.out.flagstat, by: [0])
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV_LIBRARY.out.versions.first())

    //
    // MODULE: BigWig coverage tracks
    //
    UCSC_BEDGRAPHTOBIGWIG_LIBRARY (
        BEDTOOLS_GENOMECOV_LIBRARY.out.bedgraph,
        PREPARE_GENOME.out.chrom_sizes
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG_LIBRARY.out.versions.first())

    ch_deeptoolsplotprofile_multiqc = Channel.empty()
    if (!params.skip_plot_profile) {
        //
        // MODULE: deepTools matrix generation for plotting over full transcript length
        //
        DEEPTOOLS_COMPUTEMATRIX_SCALE_REGIONS (
            UCSC_BEDGRAPHTOBIGWIG_LIBRARY.out.bigwig,
            PREPARE_GENOME.out.gene_bed
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_COMPUTEMATRIX_SCALE_REGIONS.out.versions.first())

        //
        // MODULE: deepTools matrix generation for plotting at TSS point
        //
        DEEPTOOLS_COMPUTEMATRIX_REFERENCE_POINT (
            UCSC_BEDGRAPHTOBIGWIG_LIBRARY.out.bigwig,
            PREPARE_GENOME.out.tss_bed
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_COMPUTEMATRIX_REFERENCE_POINT.out.versions.first())

        //
        // MODULE: deepTools profile plots
        //
        DEEPTOOLS_PLOTPROFILE (
            DEEPTOOLS_COMPUTEMATRIX_SCALE_REGIONS.out.matrix
        )
        ch_deeptoolsplotprofile_multiqc = DEEPTOOLS_PLOTPROFILE.out.table
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTPROFILE.out.versions.first())

        //
        // MODULE: deepTools heatmaps
        //
        DEEPTOOLS_PLOTHEATMAP (
            DEEPTOOLS_COMPUTEMATRIX_REFERENCE_POINT.out.matrix
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTHEATMAP.out.versions.first())
    }

    //
    //  Create channels: [ val(meta), [bam], [bai] ]
    //
    FILTER_BAM_BAMTOOLS
        .out
        .bam
        .join (FILTER_BAM_BAMTOOLS.out.bai, by: [0])
        .set { ch_bam_bai }

    //
    // plotFingerprint for IP and control together
    //
    ch_deeptoolsplotfingerprint_multiqc = Channel.empty()
    if (!params.skip_plot_fingerprint) {
        //
        // MODULE: deepTools plotFingerprint QC
        //
        DEEPTOOLS_PLOTFINGERPRINT (
            ch_bam_bai
        )
        ch_deeptoolsplotfingerprint_multiqc = DEEPTOOLS_PLOTFINGERPRINT.out.matrix
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTFINGERPRINT.out.versions.first())
    }

    //
    // MERGED LIBRARY PEAK ANALYSIS
    //
    ch_macs_gsize                             = params.macs_gsize
    ch_macs2_peaks_library                    = Channel.empty()
    ch_custompeaks_frip_multiqc_library       = Channel.empty()
    ch_custompeaks_count_multiqc_library      = Channel.empty()
    ch_plothomerannotatepeaks_multiqc_library = Channel.empty()
    ch_macs2_consensus_bed_library            = Channel.empty()
    ch_subreadfeaturecounts_multiqc_library   = Channel.empty()
    if (!params.macs_gsize) {
        //
        // MODULE: Estimate macs genome size
        //
        KHMER_UNIQUEKMERS (
            PREPARE_GENOME.out.fasta,
            params.read_length
        )
        ch_macs_gsize = KHMER_UNIQUEKMERS.out.kmers.map { it.text.trim() }
    }

    // Create channels: [ meta, bam, empty_list ]
    ch_bam_bai
        .map { meta, bam, bai -> [ meta , bam, [] ] }
        .set { ch_bam_library }

    //
    // MODULE: Call peaks with MACS2
    //
    MACS2_CALLPEAK_LIBRARY (
        ch_bam_library,
        ch_macs_gsize
    )
    ch_versions = ch_versions.mix(MACS2_CALLPEAK_LIBRARY.out.versions.first())

    //
    // Filter out samples with 0 MACS2 peaks called
    //
    MACS2_CALLPEAK_LIBRARY
        .out
        .peak
        .filter { meta, peaks -> peaks.size() > 0 }
        .set { ch_macs2_peaks_library }

    // Create channels: [ meta, bam, peaks ]
    ch_bam_library
        .join(ch_macs2_peaks_library, by: [0])
        .map {
            it ->
                [ it[0], it[1], it[3] ]
        }
        .set { ch_bam_peaks_library }

    //
    // MODULE: Calculate FRiP score
    //
    FRIP_SCORE_LIBRARY (
        ch_bam_peaks_library
    )
    ch_versions = ch_versions.mix(FRIP_SCORE_LIBRARY.out.versions.first())

    // Create channels: [ meta, peaks, frip ]
    ch_bam_peaks_library
        .join(FRIP_SCORE_LIBRARY.out.txt, by: [0])
        .map {
            it ->
                [ it[0], it[2], it[3] ]
        }
        .set { ch_bam_peak_frip_library }

    //
    // MODULE: FRiP score custom content for MultiQC
    //
    MULTIQC_CUSTOM_PEAKS_LIBRARY (
        ch_bam_peak_frip_library,
        ch_merged_library_peak_count_header,
        ch_merged_library_frip_score_header
    )
    ch_custompeaks_frip_multiqc_library  = MULTIQC_CUSTOM_PEAKS_LIBRARY.out.frip
    ch_custompeaks_count_multiqc_library = MULTIQC_CUSTOM_PEAKS_LIBRARY.out.count

    if (!params.skip_peak_annotation) {
        //
        // MODULE: Annotate peaks with MACS2
        //
        HOMER_ANNOTATEPEAKS_MACS2_LIBRARY (
            ch_macs2_peaks_library,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gtf
        )
        ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS_MACS2_LIBRARY.out.versions.first())

        if (!params.skip_peak_qc) {
            //
            // MODULE: MACS2 QC plots with R
            //
            PLOT_MACS2_QC_LIBRARY (
                ch_macs2_peaks_library.collect{it[1]}
            )
            ch_versions = ch_versions.mix(PLOT_MACS2_QC_LIBRARY.out.versions)

            //
            // MODULE: Peak annotation QC plots with R
            //
            PLOT_HOMER_ANNOTATEPEAKS_LIBRARY (
                HOMER_ANNOTATEPEAKS_MACS2_LIBRARY.out.txt.collect{it[1]},
                ch_merged_library_peak_annotation_header,
                "_peaks.annotatePeaks.txt"
            )
            ch_plothomerannotatepeaks_multiqc_library = PLOT_HOMER_ANNOTATEPEAKS_LIBRARY.out.tsv
            ch_versions = ch_versions.mix(PLOT_HOMER_ANNOTATEPEAKS_LIBRARY.out.versions)
        }
    }

    //
    //  Consensus peaks analysis
    //
    ch_macs2_consensus_bed_library       = Channel.empty()
    ch_deseq2_library_pca_multiqc        = Channel.empty()
    ch_deseq2_library_clustering_multiqc = Channel.empty()

    if (!params.skip_consensus_peaks) {
        // Create channels: [ meta , [ peaks ] ]
        // Where meta = [ id:consensus_peaks, multiple_groups:true/false, replicates_exist:true/false ]
        ch_macs2_peaks_library
            .map {
                meta, peak ->
                    [ '', meta.id.split('_')[0..-2].join('_'), peak ]
            }
            .groupTuple()
            .map {
                key, groups, peaks ->
                    [
                        groups.groupBy().collectEntries { [(it.key) : it.value.size()] },
                        peaks
                    ] }
            .map {
                groups, peaks ->
                    def meta_library = [:]
                    meta_library.id = 'consensus_peaks'
                    meta_library.multiple_groups = groups.size() > 1
                    meta_library.replicates_exist = groups.max { groups.value }.value > 1
                    [ meta_library, peaks ] }
            .set { ch_consensus_peaks_library }

        //
        // MODULE: Generate consensus peaks across samples
        //
        MACS2_CONSENSUS_LIBRARY (
            ch_consensus_peaks_library
        )
        ch_macs2_consensus_bed_library = MACS2_CONSENSUS_LIBRARY.out.bed
        ch_versions = ch_versions.mix(MACS2_CONSENSUS_LIBRARY.out.versions)

        if (!params.skip_peak_annotation) {
            //
            // MODULE: Annotate consensus peaks
            //
            HOMER_ANNOTATEPEAKS_CONSENSUS_LIBRARY (
                MACS2_CONSENSUS_LIBRARY.out.bed,
                PREPARE_GENOME.out.fasta,
                PREPARE_GENOME.out.gtf
            )
            ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS_CONSENSUS_LIBRARY.out.versions)
        }

        // Create channels: [ antibody, [ ip_bams ] ]
        ch_bam_library
            .map {
                meta, bam, emptylist ->
                    [ 'consensus_peaks',  meta, bam ]
            }
            .groupTuple()
            .set { ch_consensus_bams_library }

        // Create channels: [ meta, [ bams ], saf ]
        MACS2_CONSENSUS_LIBRARY
            .out
            .saf
            .map {
                meta, saf ->
                    [ meta.id, meta, saf ]
            }
            .set { ch_saf_library }

        ch_bam_library
            .map { meta, bam, emptylist -> [ 'consensus_peaks',  meta, bam ] }
            .groupTuple()
            .map {
                it ->
                    [ it[0], it[1][0], it[2].flatten().sort() ]
            }
            .join(ch_saf_library)
            .map {
                it ->
                    def fmeta_library = it[1]
                    fmeta_library['id'] = it[3]['id']
                    fmeta_library['replicates_exist'] = it[3]['replicates_exist']
                    fmeta_library['multiple_groups']  = it[3]['multiple_groups']
                    [ fmeta_library, it[2], it[4] ] }
            .set { ch_saf_bams_library }

        SUBREAD_FEATURECOUNTS_LIBRARY (
            ch_saf_bams_library
        )
        ch_subreadfeaturecounts_multiqc_library = SUBREAD_FEATURECOUNTS_LIBRARY.out.summary
        ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS_LIBRARY.out.versions.first())

        if (!params.skip_deseq2_qc) {
            //
            // MODULE: Generate QC plots with DESeq2
            //
            DESEQ2_QC_LIBRARY (
                SUBREAD_FEATURECOUNTS_LIBRARY.out.counts,
                ch_merged_library_deseq2_pca_header,
                ch_merged_library_deseq2_clustering_header
            )
            ch_deseq2_library_pca_multiqc        = DESEQ2_QC_LIBRARY.out.pca_multiqc
            ch_deseq2_library_clustering_multiqc = DESEQ2_QC_LIBRARY.out.dists_multiqc
        }
    }

    // Create channel: [ val(meta), bam, bai, peak_file ]
    MARK_DUPLICATES_PICARD_LIBRARY
        .out
        .bam
        .join(MARK_DUPLICATES_PICARD_LIBRARY.out.bai, by: [0])
        .join(ch_macs2_peaks_library, by: [0])
        .set { ch_bam_peak }

    if (!params.skip_ataqv) {
        //
        // MODULE: ATAQV QC
        //
        ATAQV_ATAQV (
            ch_bam_peak,
            'NA',
            params.mito_name,
            PREPARE_GENOME.out.tss_bed,
            [],
            PREPARE_GENOME.out.autosomes
        )
        ch_versions = ch_versions.mix(ATAQV_ATAQV.out.versions)

        ATAQV_MKARV (
            ATAQV_ATAQV.out.json.collect{it[1]}
        )
        ch_versions = ch_versions.mix(ATAQV_MKARV.out.versions)
    }

    //
    // MERGED REPLICATES ANALYSIS
    //
    ch_ucsc_bedgraphtobigwig_replicate_bigwig    = Channel.empty()
    ch_macs2_peaks_replicate                     = Channel.empty()
    ch_macs2_consensus_bed_replicate             = Channel.empty()
    ch_mark_duplicates_picard_replicate_stats    = Channel.empty()
    ch_mark_duplicates_picard_replicate_flagstat = Channel.empty()
    ch_mark_duplicates_picard_replicate_idxstats = Channel.empty()
    ch_mark_duplicates_picard_replicate_metrics  = Channel.empty()
    ch_custompeaks_frip_multiqc_replicate        = Channel.empty()
    ch_custompeaks_count_multiqc_replicate       = Channel.empty()
    ch_plothomerannotatepeaks_multiqc_replicate  = Channel.empty()
    ch_subreadfeaturecounts_multiqc_replicate    = Channel.empty()
    ch_deseq2_replicate_pca_multiqc              = Channel.empty()
    ch_deseq2_replicate_clustering_multiqc       = Channel.empty()
    if (!params.skip_merge_replicates) {
        //
        // MERGE REPLICATE BAM
        // MODULE: Merge replicate BAM files
        //
        FILTER_BAM_BAMTOOLS
            .out
            .bam
            .map {
                meta, bam ->
                    def fmeta = meta.findAll { it.key != 'read_group' }
                    fmeta.id = fmeta.id.split('_')[0..-2].join('_')
                    [ fmeta, fmeta.id, bam ]
            }
            .groupTuple(by: [0])
            .map {
                meta, ids, bams ->
                    meta.replicates_exist = ids.size() > 1
                    meta.multiple_bams    = bams.size() > 1
                    [ meta, bams.flatten() ]
            }
            .set { ch_merged_library_rm_orphan_bam_merged_replicate }

        PICARD_MERGESAMFILES_REPLICATE (
            ch_merged_library_rm_orphan_bam_merged_replicate
        )
        ch_versions = ch_versions.mix(PICARD_MERGESAMFILES_REPLICATE.out.versions.first().ifEmpty(null))

        //
        // SUBWORKFLOW: Mark duplicates & filter BAM files after merging
        //
        MARK_DUPLICATES_PICARD_REPLICATE (
            PICARD_MERGESAMFILES_REPLICATE.out.bam
        )
        ch_mark_duplicates_picard_replicate_stats    = MARK_DUPLICATES_PICARD_REPLICATE.out.stats
        ch_mark_duplicates_picard_replicate_flagstat = MARK_DUPLICATES_PICARD_REPLICATE.out.flagstat
        ch_mark_duplicates_picard_replicate_idxstats = MARK_DUPLICATES_PICARD_REPLICATE.out.idxstats
        ch_mark_duplicates_picard_replicate_metrics  = MARK_DUPLICATES_PICARD_REPLICATE.out.metrics

        //
        // MODULE: BigWig coverage tracks
        //
        BEDTOOLS_GENOMECOV_REPLICATE (
            MARK_DUPLICATES_PICARD_REPLICATE.out.bam.join(MARK_DUPLICATES_PICARD_REPLICATE.out.flagstat, by: [0])
        )
        ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV_REPLICATE.out.versions.first())

        //
        // MODULE: BigWig coverage tracks
        //
        UCSC_BEDGRAPHTOBIGWIG_REPLICATE (
            BEDTOOLS_GENOMECOV_REPLICATE.out.bedgraph,
            PREPARE_GENOME.out.chrom_sizes
        )
        ch_ucsc_bedgraphtobigwig_replicate_bigwig = UCSC_BEDGRAPHTOBIGWIG_REPLICATE.out.bigwig
        ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG_REPLICATE.out.versions.first())

        //
        // MERGE REPLICATES PEAK ANALYSIS
        // Create channel: [ val(meta), bam, bai, peak_file ]
        MARK_DUPLICATES_PICARD_REPLICATE
            .out
            .bam
            .map {
                meta, bam ->
                    [ meta , bam, [] ]
            }
            .set { ch_bam_replicate }

        //
        // MODULE: Call peaks with MACS2
        //
        MACS2_CALLPEAK_REPLICATE (
            ch_bam_replicate,
            ch_macs_gsize
        )
        ch_versions = ch_versions.mix(MACS2_CALLPEAK_REPLICATE.out.versions.first())

        //
        // Filter out samples with 0 MACS2 peaks called
        //
        MACS2_CALLPEAK_REPLICATE
            .out
            .peak
            .filter { meta, peaks -> peaks.size() > 0 }
            .set { ch_macs2_peaks_replicate }

        // Create channels: [ meta, bam, peaks ]
        ch_bam_replicate
            .join(ch_macs2_peaks_replicate, by: [0])
            .map {
                it ->
                    [ it[0], it[1], it[3] ]
            }
            .set { ch_bam_peaks_replicate }

        //
        // MODULE: Calculate FRiP score
        //
        FRIP_SCORE_REPLICATE (
            ch_bam_peaks_replicate
        )
        ch_versions = ch_versions.mix(FRIP_SCORE_REPLICATE.out.versions.first())

        // Create channels: [ meta, peaks, frip ]
        ch_bam_peaks_replicate
            .join(FRIP_SCORE_REPLICATE.out.txt, by: [0])
            .map {
                it ->
                    [ it[0], it[2], it[3] ]
            }
            .set { ch_bam_peak_frip_replicate }

        //
        // MODULE: FRiP score custom content for MultiQC
        //
        MULTIQC_CUSTOM_PEAKS_REPLICATE (
            ch_bam_peak_frip_replicate,
            ch_merged_replicate_peak_count_header,
            ch_merged_replicate_frip_score_header
        )
        ch_custompeaks_frip_multiqc_replicate  = MULTIQC_CUSTOM_PEAKS_REPLICATE.out.frip
        ch_custompeaks_count_multiqc_replicate = MULTIQC_CUSTOM_PEAKS_REPLICATE.out.count

        if (!params.skip_peak_annotation) {
            //
            // MODULE: Annotate peaks with MACS2
            //
            HOMER_ANNOTATEPEAKS_MACS2_REPLICATE (
                ch_macs2_peaks_replicate,
                PREPARE_GENOME.out.fasta,
                PREPARE_GENOME.out.gtf
            )
            ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS_MACS2_REPLICATE.out.versions.first())

            if (!params.skip_peak_qc) {
                //
                // MODULE: MACS2 QC plots with R
                //
                PLOT_MACS2_QC_REPLICATE (
                    ch_macs2_peaks_replicate.collect{it[1]}
                )
                ch_versions = ch_versions.mix(PLOT_MACS2_QC_REPLICATE.out.versions)

                //
                // MODULE: Peak annotation QC plots with R
                //
                PLOT_HOMER_ANNOTATEPEAKS_REPLICATE (
                    HOMER_ANNOTATEPEAKS_MACS2_REPLICATE.out.txt.collect{it[1]},
                    ch_merged_replicate_peak_annotation_header,
                    "_peaks.annotatePeaks.txt"
                )
                ch_plothomerannotatepeaks_multiqc_replicate = PLOT_HOMER_ANNOTATEPEAKS_REPLICATE.out.tsv
                ch_versions = ch_versions.mix(PLOT_HOMER_ANNOTATEPEAKS_REPLICATE.out.versions)
            }
        }

        //
        //  Consensus peaks analysis
        //
        if (!params.skip_consensus_peaks) {
            // Create channel: [ meta , [ peaks ] ]
            // Where meta = [ id:consensus_peaks, multiple_groups:true/false, replicates_exist:true/false ]
            ch_macs2_peaks_replicate
                .map {
                    meta, peak ->
                        [ '', meta.id, meta, peak ]
                }
                .groupTuple()
                .map {
                    key, groups, metas, peaks ->
                        def meta_replicate = [:]
                        meta_replicate.id = 'consensus_peaks'
                        meta_replicate.replicates_exist = metas[0].replicates_exist
                        meta_replicate.multiple_groups = groups.groupBy().collectEntries { [(it.key) : it.value.size()] }.size() > 1
                        [ meta_replicate, peaks ]
                }
                .set { ch_consensus_peaks_replicate }

            //
            // MODULE: Generate consensus peaks across samples
            //
            MACS2_CONSENSUS_REPLICATE (
                ch_consensus_peaks_replicate
            )
            ch_macs2_consensus_bed_replicate = MACS2_CONSENSUS_REPLICATE.out.bed
            ch_versions = ch_versions.mix(MACS2_CONSENSUS_REPLICATE.out.versions)

            if (!params.skip_peak_annotation) {
                //
                // MODULE: Annotate consensus peaks
                //
                HOMER_ANNOTATEPEAKS_CONSENSUS_REPLICATE (
                    MACS2_CONSENSUS_REPLICATE.out.bed,
                    PREPARE_GENOME.out.fasta,
                    PREPARE_GENOME.out.gtf
                )
                ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS_CONSENSUS_REPLICATE.out.versions)
            }

            // Create channels: [ antibody, [ ip_bams ] ]
            ch_bam_replicate
                .map {
                    meta, bam, emptylist ->
                        [ 'consensus_peaks',  meta, bam ]
                }
                .groupTuple()
                .set { ch_consensus_bams_replicate }

            // Create channels: [ meta, [ bams ], saf ]
            MACS2_CONSENSUS_REPLICATE
                .out
                .saf
                .map {
                    meta, saf ->
                        [ meta.id, meta, saf ]
                }
                .set { ch_saf_replicate }

            ch_bam_replicate
                .map {
                    meta, bam, emptylist ->
                        [ 'consensus_peaks',  meta, bam ]
                }
                .groupTuple()
                .map {
                    it ->
                        [ it[0], it[1][0], it[2].flatten().sort() ]
                }
                .join(ch_saf_replicate)
                .map {
                    it ->
                        def fmeta_replicate = it[1]
                        fmeta_replicate['id'] = it[3]['id']
                        fmeta_replicate['replicates_exist'] = it[3]['replicates_exist']
                        fmeta_replicate['multiple_groups']  = it[3]['multiple_groups']
                        [ fmeta_replicate, it[2], it[4] ] }
                .set { ch_saf_bams_replicate }

            //
            // MODULE: Quantify peaks across samples with featureCounts
            //
            SUBREAD_FEATURECOUNTS_REPLICATE (
                ch_saf_bams_replicate
            )
            ch_subreadfeaturecounts_multiqc_replicate = SUBREAD_FEATURECOUNTS_REPLICATE.out.summary
            ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS_REPLICATE.out.versions.first())

            if (!params.skip_deseq2_qc) {
                //
                // MODULE: Generate QC plots with DESeq2
                //
                DESEQ2_QC_REPLICATE (
                    SUBREAD_FEATURECOUNTS_REPLICATE.out.counts,
                    ch_merged_replicate_deseq2_pca_header,
                    ch_merged_replicate_deseq2_clustering_header
                )
                ch_deseq2_replicate_pca_multiqc        = DESEQ2_QC_REPLICATE.out.pca_multiqc
                ch_deseq2_replicate_clustering_multiqc = DESEQ2_QC_REPLICATE.out.dists_multiqc
            }
        }
    }

    //
    // MODULE: Create IGV session
    //
    if (!params.skip_igv) {
        IGV (
            PREPARE_GENOME.out.fasta,
            UCSC_BEDGRAPHTOBIGWIG_LIBRARY.out.bigwig.collect{it[1]}.ifEmpty([]),
            ch_macs2_peaks_library.collect{it[1]}.ifEmpty([]),
            ch_macs2_consensus_bed_library.collect{it[1]}.ifEmpty([]),
            ch_ucsc_bedgraphtobigwig_replicate_bigwig.collect{it[1]}.ifEmpty([]),
            ch_macs2_peaks_replicate.collect{it[1]}.ifEmpty([]),
            ch_macs2_consensus_bed_replicate.collect{it[1]}.ifEmpty([]),
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
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowAtacseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        methods_description    = WorkflowAtacseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
        ch_methods_description = Channel.value(methods_description)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),

            FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
            FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]),
            FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]),

            ch_samtools_stats.collect{it[1]}.ifEmpty([]),
            ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
            ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),

            MARK_DUPLICATES_PICARD_LIBRARY.out.stats.collect{it[1]}.ifEmpty([]),
            MARK_DUPLICATES_PICARD_LIBRARY.out.flagstat.collect{it[1]}.ifEmpty([]),
            MARK_DUPLICATES_PICARD_LIBRARY.out.idxstats.collect{it[1]}.ifEmpty([]),
            MARK_DUPLICATES_PICARD_LIBRARY.out.metrics.collect{it[1]}.ifEmpty([]),

            FILTER_BAM_BAMTOOLS.out.stats.collect{it[1]}.ifEmpty([]),
            FILTER_BAM_BAMTOOLS.out.flagstat.collect{it[1]}.ifEmpty([]),
            FILTER_BAM_BAMTOOLS.out.idxstats.collect{it[1]}.ifEmpty([]),
            ch_picardcollectmultiplemetrics_multiqc.collect{it[1]}.ifEmpty([]),

            ch_preseq_multiqc.collect{it[1]}.ifEmpty([]),

            ch_deeptoolsplotprofile_multiqc.collect{it[1]}.ifEmpty([]),
            ch_deeptoolsplotfingerprint_multiqc.collect{it[1]}.ifEmpty([]),

            ch_custompeaks_frip_multiqc_library.collect{it[1]}.ifEmpty([]),
            ch_custompeaks_count_multiqc_library.collect{it[1]}.ifEmpty([]),
            ch_plothomerannotatepeaks_multiqc_library.collect().ifEmpty([]),
            ch_subreadfeaturecounts_multiqc_library.collect{it[1]}.ifEmpty([]),

            ch_mark_duplicates_picard_replicate_stats.collect{it[1]}.ifEmpty([]),
            ch_mark_duplicates_picard_replicate_flagstat.collect{it[1]}.ifEmpty([]),
            ch_mark_duplicates_picard_replicate_idxstats.collect{it[1]}.ifEmpty([]),
            ch_mark_duplicates_picard_replicate_metrics.collect{it[1]}.ifEmpty([]),

            ch_custompeaks_frip_multiqc_replicate.collect{it[1]}.ifEmpty([]),
            ch_custompeaks_count_multiqc_replicate.collect{it[1]}.ifEmpty([]),
            ch_plothomerannotatepeaks_multiqc_replicate.collect().ifEmpty([]),
            ch_subreadfeaturecounts_multiqc_replicate.collect{it[1]}.ifEmpty([]),

            ch_deseq2_library_pca_multiqc.collect().ifEmpty([]),
            ch_deseq2_library_clustering_multiqc.collect().ifEmpty([]),
            ch_deseq2_replicate_pca_multiqc.collect().ifEmpty([]),
            ch_deseq2_replicate_clustering_multiqc.collect().ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }

    if (params.hook_url) {
        NfcoreTemplate.adaptivecard(workflow, params, summary_params, projectDir, log)
    }

    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
