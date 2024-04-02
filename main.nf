#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/atacseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/atacseq
    Website: https://nf-co.re/atacseq
    Slack  : https://nfcore.slack.com/channels/atacseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ATACSEQ                 } from './workflows/atacseq'
include { PREPARE_GENOME          } from './subworkflows/local/prepare_genome'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_atacseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_atacseq_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_atacseq_pipeline'
include { getMacsGsize            } from './subworkflows/local/utils_nfcore_atacseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta         = getGenomeAttribute('fasta')
params.bwa_index     = getGenomeAttribute('bwa')
params.bowtie2_index = getGenomeAttribute('bowtie2')
params.chromap_index = getGenomeAttribute('chromap')
params.star_index    = getGenomeAttribute('star')
params.gtf           = getGenomeAttribute('gtf')
params.gff           = getGenomeAttribute('gff')
params.gene_bed      = getGenomeAttribute('gene_bed')
params.tss_bed       = getGenomeAttribute('tss_bed')
params.blacklist     = getGenomeAttribute('blacklist')
params.mito_name     = getGenomeAttribute('mito_name')
params.macs_gsize    = getMacsGsize(params)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline
//
workflow NFCORE_ATACSEQ {

    main:
    ch_versions = Channel.empty()

    // SUBWORKFLOW: Prepare genome files
    PREPARE_GENOME (
        params.genome,
        params.genomes,
        params.aligner,
        params.fasta,
        params.gtf,
        params.gff,
        params.blacklist,
        params.gene_bed,
        params.tss_bed,
        params.mito_name,
        params.keep_mito,
        params.bwa_index,
        params.bowtie2_index,
        params.chromap_index,
        params.star_index,
        params.macs_gsize,
        params.read_length
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // WORKFLOW: Run nf-core/atacseq workflow
    //
    ch_samplesheet = Channel.value(file(params.input, checkIfExists: true))

    ATACSEQ (
        params.input,
        ch_versions,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fai,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.gene_bed,
        PREPARE_GENOME.out.tss_bed,
        PREPARE_GENOME.out.chrom_sizes,
        PREPARE_GENOME.out.filtered_bed,
        PREPARE_GENOME.out.bwa_index,
        PREPARE_GENOME.out.bowtie2_index,
        PREPARE_GENOME.out.chromap_index,
        PREPARE_GENOME.out.star_index,
        PREPARE_GENOME.out.autosomes,
        PREPARE_GENOME.out.macs_gsize
    )

    emit:
    multiqc_report = ATACSEQ.out.multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions                // channel: [version1, version2, ...]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_ATACSEQ ()

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_ATACSEQ.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
