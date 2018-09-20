#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/atacseq
========================================================================================
 nf-core/atacseq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/atacseq
----------------------------------------------------------------------------------------
*/

//nf-core/atacseq v${manifest.pipelineVersion}
def helpMessage() {
    log.info"""
    =========================================

    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/atacseq --design design.csv --genome GRCh37 -profile standard,docker

    Mandatory arguments:
      --design                      Comma-separted file containing information about the samples in the experiment (see README.md)
      --genome                      Name of iGenomes reference
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: standard, conda, docker, singularity, awsbatch, test

    Options:

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference
      --bwa_index                   Path to BWA index
      --gtf                         Path to GTF file (Ensembl format)
      --bed12                       Path to bed12 file
      --mito_name                   Name of Mitochondrial chomosome in genome fasta (e.g. chrM). Reads aligning to this contig are filtered out
      --macs_gsize                  Effective genome size parameter required by MACS2. It only works when --genome is set as GRCh37, GRCm38, BDGP6 and WBcel235
      --blacklist                   Path to blacklist regions (.BED format), used for filtering out called peaks. Note that --blacklist_filtering is required
      --saveReference               Save the generated reference files in the Results directory.

    Alignments
      --blacklist_filtering         Filter ENCODE blacklisted regions from ATAC-seq peaks. It only works when --genome is set as GRCh37 or GRCm38
      --allow_multi_align           Reads mapping to multiple places are not filtered from alignments
      --keep_duplicates             Duplicate reads are not filtered from alignments
      --saveAlignedIntermediates    Save the intermediate BAM files from the Alignment step  - not done by default

    Trimming
      --adapter                     3' adapter sequence trimmed by cutadapt. Default: 'CTGTCTCTTATA'
      --notrim                      Specifying --notrim will skip the adapter trimming step.

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Configurable variables
params.name = false
params.genome = false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.bed12 = params.genome ? params.genomes[ params.genome ].bed12 ?: false : false
params.mito_name = params.genome ? params.genomes[ params.genome ].mito_name ?: false : false
params.macs_gsize = params.genome ? params.genomes[ params.genome ].macs_gsize ?: false : false
params.blacklist = params.genome ? params.genomes[ params.genome ].blacklist ?: false : false
params.adapter = false
params.multiqc_config = "$baseDir/assets/multiqc_config.yaml"
params.bamtools_filter_config = "$baseDir/conf/bamtools_filter_pe.json"
params.email = false
params.plaintext_email = false

multiqc_config = file(params.multiqc_config)
bamtools_filter_config = file(params.bamtools_filter_config)
output_docs = file("$baseDir/docs/output.md")
wherearemyfiles = file("$baseDir/assets/where_are_my_files.txt")

// Preset adapter trimming options
adapter = 'CTGTCTCTTATA'
if (params.adapter){
    adapter = params.adapter
}

// Check if macs_gsize is valid
def RUN_macs = false
if (params.macs_gsize != ""){ RUN_macs = true }

// Validate inputs
if( params.bwa_index ){
    bwa_index = Channel
        .fromPath(params.bwa_index)
        .ifEmpty { exit 1, "BWA index not found: ${params.bwa_index}" }
} else if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
} else {
    exit 1, "No reference genome specified!"
}

if( params.gtf ){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .into { gtf_1; gtf_2 }
}

if( params.bed12 ){
    bed12 = Channel
        .fromPath(params.bed12)
        .ifEmpty { exit 1, "BED12 annotation file not found: ${params.bed12}" }
        .into { bed_1; bed_2 }
}

if ( params.blacklist_filtering ){
    blacklist = file(params.blacklist)
    if( !blacklist.exists() ) exit 1, "Blacklist file not found: ${params.blacklist}"
}

// AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Check workDir/outdir paths to be S3 buckets if running on AWSBatch
// related: https://github.com/nextflow-io/nextflow/issues/813
if( workflow.profile == 'awsbatch') {
    if(!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

/*
 * Create a channel for input read files
 */
if( params.design ){
   Channel
     .fromPath(params.design)
     .ifEmpty { exit 1, "Design file not found: ${params.design}" }
     .splitCsv(header:true, sep:',')
     .map { row -> [ [row.sample,"R"+row.replicate,"L"+row.run].join("_"),
                      row.sample,
                      row.replicate,
                      row.run,
                      file(row.fastq_1),
                      file(row.fastq_2) ] }
     .into { design_replicates_exist_ch;
             design_multiple_samples_ch;
             design_raw_fastqc_ch;
             design_raw_fastqscreen_ch;
             design_cutadapt_ch }
}

// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'


======================================================="""
//nf-core/atacseq v${manifest.pipelineVersion}"
def summary = [:]
summary['Pipeline Name']          = 'nf-core/atacseq'
summary['Pipeline Version']       = manifest.pipelineVersion
summary['Run Name'] = custom_runName ?: workflow.runName
summary['Design File']            = params.design
summary['Genome']                 = params.genome
if(params.bwa_index)  summary['BWA Index'] = params.bwa_index
else if(params.fasta) summary['Fasta Ref'] = params.fasta
summary['GTF File']               = params.gtf
summary['BED12 File']             = params.bed12
summary['Mitochondrial Contig']   = params.mito_name
summary['MACS2 Genome Size']      = params.macs_gsize
if( params.blacklist_filtering ) summary['Blacklist BED'] = params.blacklist
summary['Save Reference']         = params.saveReference ? 'Yes' : 'No'
summary['Blacklist Filtering']    = params.blacklist_filtering ? 'Yes' : 'No'
summary['Allow Multi-mapped']     = params.allow_multi_align ? 'Yes' : 'No'
summary['Keep Duplicates']        = params.keep_duplicates ? 'Yes' : 'No'
summary['Save Aligned Intermeds'] = params.saveAlignedIntermediates ? 'Yes' : 'No'
summary['Trimming Step']          = params.notrim ? 'Yes' : 'No'
summary['Adapter']                = adapter
summary['Max Memory']             = params.max_memory
summary['Max CPUs']               = params.max_cpus
summary['Max Time']               = params.max_time
summary['Output dir']             = params.outdir
summary['Working dir']            = workflow.workDir
summary['Container Engine']       = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']           = "$HOME"
summary['Current user']           = "$USER"
summary['Current path']           = "$PWD"
summary['Working dir']            = workflow.workDir
summary['Output dir']             = params.outdir
summary['Script dir']             = workflow.projectDir
summary['Config Profile']         = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']          = params.awsregion
   summary['AWS Queue']           = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(23)}: $v" }.join("\n")
log.info "========================================="

// // Show a big warning message if we're not running MACS
// if (!RUN_macs){
//     def warnstring = params.genome ? "Reference '${params.genome}' not supported by" : 'No reference supplied for'
//     log.warn "=======================================================\n" +
//              "  WARNING! $warnstring MACS, ngs_plot\n" +
//              "  and annotation. Steps for MACS, ngs_plot and annotation\n" +
//              "  will be skipped. Use '--genome GRCh37' or '--genome GRCm38'\n" +
//              "  to run these steps.\n" +
//              "==============================================================="
// }

// def create_workflow_summary(summary) {
//
//     def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
//     yaml_file.text  = """
//     id: 'nf-core-atacseq-summary'
//     description: " - this information is collected when the pipeline is started."
//     section_name: 'nf-core/atacseq Workflow Summary'
//     section_href: 'https://github.com/nf-core/atacseq'
//     plot_type: 'html'
//     data: |
//         <dl class=\"dl-horizontal\">
// ${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
//         </dl>
//     """.stripIndent()
//
//    return yaml_file
// }
//
// /*
//  * Parse software version numbers
//  */
// process get_software_versions {
//
//     output:
//     file 'software_versions_mqc.yaml' into software_versions_yaml
//
//     script:
//     """
//     echo $manifest.pipelineVersion > v_pipeline.txt
//     echo $workflow.nextflow.version > v_nextflow.txt
//     fastqc --version > v_fastqc.txt
//     multiqc --version > v_multiqc.txt
//     scrape_software_versions.py > software_versions_mqc.yaml
//     """
// }
//
//
// /*
//  * PREPROCESSING - Build BWA index
//  */
// if(!params.bwa_index && fasta){
//     process makeBWAindex {
//         tag fasta
//         publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
//                    saveAs: { params.saveReference ? it : null }, mode: 'copy'
//
//         input:
//         file fasta from fasta
//
//         output:
//         file "BWAIndex" into bwa_index
//
//         script:
//         """
//         bwa index -a bwtsw $fasta
//         mkdir BWAIndex && mv ${fasta}* BWAIndex
//         """
//     }
// }
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
// /*
//  * STEP 1 - FastQC
//  */
// process fastqc {
//     tag "$name"
//     publishDir "${params.outdir}/fastqc", mode: 'copy',
//         saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}
//
//     input:
//     set val(name), file(reads) from read_files_fastqc
//
//     output:
//     file "*_fastqc.{zip,html}" into fastqc_results
//
//     script:
//     """
//     fastqc -q $reads
//     """
// }
//
//
//
// /*
//  * STEP 2 - MultiQC
//  */
// process multiqc {
//     publishDir "${params.outdir}/MultiQC", mode: 'copy'
//
//     input:
//     file multiqc_config
//     file ('fastqc/*') from fastqc_results.collect()
//     file ('software_versions/*') from software_versions_yaml
//     file workflow_summary from create_workflow_summary(summary)
//
//     output:
//     file "*multiqc_report.html" into multiqc_report
//     file "*_data"
//
//     script:
//     rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
//     rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
//     """
//     multiqc -f $rtitle $rfilename --config $multiqc_config .
//     """
// }
//
//
//
// /*
//  * STEP 3 - Output Description HTML
//  */
// process output_documentation {
//     tag "$prefix"
//     publishDir "${params.outdir}/Documentation", mode: 'copy'
//
//     input:
//     file output_docs
//
//     output:
//     file "results_description.html"
//
//     script:
//     """
//     markdown_to_html.r $output_docs results_description.html
//     """
// }
//
//
//
// /*
//  * Completion e-mail notification
//  */
// workflow.onComplete {
//
//     // Set up the e-mail variables
//     def subject = "[nf-core/atacseq] Successful: $workflow.runName"
//     if(!workflow.success){
//       subject = "[nf-core/atacseq] FAILED: $workflow.runName"
//     }
//     def email_fields = [:]
//     email_fields['version'] = manifest.pipelineVersion
//     email_fields['runName'] = custom_runName ?: workflow.runName
//     email_fields['success'] = workflow.success
//     email_fields['dateComplete'] = workflow.complete
//     email_fields['duration'] = workflow.duration
//     email_fields['exitStatus'] = workflow.exitStatus
//     email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
//     email_fields['errorReport'] = (workflow.errorReport ?: 'None')
//     email_fields['commandLine'] = workflow.commandLine
//     email_fields['projectDir'] = workflow.projectDir
//     email_fields['summary'] = summary
//     email_fields['summary']['Date Started'] = workflow.start
//     email_fields['summary']['Date Completed'] = workflow.complete
//     email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
//     email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
//     if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
//     if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
//     if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
//     email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
//     email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
//     email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp
//
//     // Render the TXT template
//     def engine = new groovy.text.GStringTemplateEngine()
//     def tf = new File("$baseDir/assets/email_template.txt")
//     def txt_template = engine.createTemplate(tf).make(email_fields)
//     def email_txt = txt_template.toString()
//
//     // Render the HTML template
//     def hf = new File("$baseDir/assets/email_template.html")
//     def html_template = engine.createTemplate(hf).make(email_fields)
//     def email_html = html_template.toString()
//
//     // Render the sendmail template
//     def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
//     def sf = new File("$baseDir/assets/sendmail_template.txt")
//     def sendmail_template = engine.createTemplate(sf).make(smail_fields)
//     def sendmail_html = sendmail_template.toString()
//
//     // Send the HTML e-mail
//     if (params.email) {
//         try {
//           if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
//           // Try to send HTML e-mail using sendmail
//           [ 'sendmail', '-t' ].execute() << sendmail_html
//           log.info "[nf-core/atacseq] Sent summary e-mail to $params.email (sendmail)"
//         } catch (all) {
//           // Catch failures and try with plaintext
//           [ 'mail', '-s', subject, params.email ].execute() << email_txt
//           log.info "[nf-core/atacseq] Sent summary e-mail to $params.email (mail)"
//         }
//     }
//
//     // Write summary e-mail HTML to a file
//     def output_d = new File( "${params.outdir}/Documentation/" )
//     if( !output_d.exists() ) {
//       output_d.mkdirs()
//     }
//     def output_hf = new File( output_d, "pipeline_report.html" )
//     output_hf.withWriter { w -> w << email_html }
//     def output_tf = new File( output_d, "pipeline_report.txt" )
//     output_tf.withWriter { w -> w << email_txt }
//
//     log.info "[nf-core/atacseq] Pipeline Complete"
//
// }
