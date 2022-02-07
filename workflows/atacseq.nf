/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowAtacseq.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input, params.multiqc_config,
    params.fasta,
    params.gtf, params.gff, params.gene_bed,
    params.bwa_index,
    params.blacklist,
    params.bamtools_filter_pe_config, params.bamtools_filter_se_config
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Save AWS IGenomes file containing annotation version
def anno_readme = params.genomes[ params.genome ]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()


// // JSON files required by BAMTools for alignment filtering
// ch_bamtools_filter_se_config = file(params.bamtools_filter_se_config, checkIfExists: true)
// ch_bamtools_filter_pe_config = file(params.bamtools_filter_pe_config, checkIfExists: true)

// // Header files for MultiQC
// ch_mlib_peak_count_header = file("$baseDir/assets/multiqc/mlib_peak_count_header.txt", checkIfExists: true)
// ch_mlib_frip_score_header = file("$baseDir/assets/multiqc/mlib_frip_score_header.txt", checkIfExists: true)
// ch_mlib_peak_annotation_header = file("$baseDir/assets/multiqc/mlib_peak_annotation_header.txt", checkIfExists: true)
// ch_mlib_deseq2_pca_header = file("$baseDir/assets/multiqc/mlib_deseq2_pca_header.txt", checkIfExists: true)
// ch_mlib_deseq2_clustering_header = file("$baseDir/assets/multiqc/mlib_deseq2_clustering_header.txt", checkIfExists: true)

// ch_mrep_peak_count_header = file("$baseDir/assets/multiqc/mrep_peak_count_header.txt", checkIfExists: true)
// ch_mrep_frip_score_header = file("$baseDir/assets/multiqc/mrep_frip_score_header.txt", checkIfExists: true)
// ch_mrep_peak_annotation_header = file("$baseDir/assets/multiqc/mrep_peak_annotation_header.txt", checkIfExists: true)
// ch_mrep_deseq2_pca_header = file("$baseDir/assets/multiqc/mrep_deseq2_pca_header.txt", checkIfExists: true)
// ch_mrep_deseq2_clustering_header = file("$baseDir/assets/multiqc/mrep_deseq2_clustering_header.txt", checkIfExists: true)

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// include { BEDTOOLS_GENOMECOV                  } from '../modules/local/bedtools_genomecov'
include { MULTIQC                             } from '../modules/local/multiqc'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK         } from '../subworkflows/local/input_check'
include { PREPARE_GENOME      } from '../subworkflows/local/prepare_genome'
// include { FILTER_BAM_BAMTOOLS } from '../subworkflows/local/filter_bam_bamtools'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

// /*
// ========================================================================================
//     RUN MAIN WORKFLOW
// ========================================================================================
// */

// // Info required for completion email and summary
// def multiqc_report = []

workflow ATACSEQ {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    PREPARE_GENOME ()
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input),
        params.seq_center
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

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

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
        )
        multiqc_report       = MULTIQC.out.report.toList()
    }
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
