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

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                             PARAMETERS                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

def helpMessage() {
    log.info"""
    =========================================
     nf-core/atacseq v${workflow.manifest.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/atacseq --design design.csv --genome GRCh37 -profile standard,docker

    Mandatory arguments:
      --design                      Comma-separted file containing information about the samples in the experiment (see docs/usage.md)
      --fasta                       Path to Fasta reference. Not mandatory when using reference in iGenomes config via --genome
      --gtf                         Path to GTF file in Ensembl format. Not mandatory when using reference in iGenomes config via --genome
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: standard, conda, docker, singularity, awsbatch, test

    Generic
      --genome                      Name of iGenomes reference
      --singleEnd                   Specifies that the input is single-end reads
      --narrowPeak                  Run MACS in narrowPeak mode. Default: broadPeak
      --fragment_size [int]         Estimated fragment size used to extend single-end reads. Default: 0

    References                      If not specified in the configuration file or you wish to overwrite any of the references
      --bwa_index                   Path to BWA index
      --bed12                       Path to bed12 file
      --mito_name                   Name of Mitochondrial chomosome in genome fasta (e.g. chrM). Reads aligning to this contig are filtered out

      --macs_gsize                  Effective genome size parameter required by MACS2. If using igenomes config, values have only been provided when --genome is set as GRCh37, GRCm38, hg19, mm10, BDGP6 and WBcel235
      --blacklist                   Path to blacklist regions (.BED format), used for filtering alignments
      --saveReference               Save the generated reference files in the Results directory

    Trimming
      --clip_r1 [int]               Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads)
      --clip_r2 [int]               Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only)
      --three_prime_clip_r1 [int]   Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed
      --three_prime_clip_r2 [int]   Instructs Trim Galore to re move bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed
      --skipTrimming                Skip the adapter trimming step
      --saveTrimmed                 Save the trimmed FastQ files in the the Results directory

    Alignments
      --keepMito                    Reads mapping to mitochondrial contig are not filtered from alignments
      --keepDups                    Duplicate reads are not filtered from alignments
      --keepMultiMap                Reads mapping to multiple locations are not filtered from alignments
      --skipMergeBySample           Do not perform alignment merging and downstream analysis at the sample-level i.e. only do this at the replicate-level
      --saveAlignedIntermediates    Save the intermediate BAM files from the alignment step - not done by default

    Other
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                SET UP CONFIGURATION VARIABLES                       -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

////////////////////////////////////////////////////
/* --         DEFAULT PARAMETER VALUES         -- */
////////////////////////////////////////////////////

// Configurable variables
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.bed12 = params.genome ? params.genomes[ params.genome ].bed12 ?: false : false
params.mito_name = params.genome ? params.genomes[ params.genome ].mito_name ?: false : false
params.macs_gsize = params.genome ? params.genomes[ params.genome ].macs_gsize ?: false : false
params.blacklist = params.genome ? params.genomes[ params.genome ].blacklist ?: false : false

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

multiqc_config_ch = Channel.fromPath(params.multiqc_config)
bamtools_filter_pe_config_ch = Channel.fromPath(params.bamtools_filter_pe_config)
bamtools_filter_se_config_ch = Channel.fromPath(params.bamtools_filter_se_config)
output_docs_ch = Channel.fromPath("$baseDir/docs/output.md")

// Header files for MultiQC custom-content
replicate_peak_count_header_ch = Channel.fromPath("$baseDir/assets/multiqc/replicate_peak_count_header.txt")
replicate_frip_score_header_ch = Channel.fromPath("$baseDir/assets/multiqc/replicate_frip_score_header.txt")
replicate_peak_annotation_header_ch = Channel.fromPath("$baseDir/assets/multiqc/replicate_peak_annotation_header.txt")
replicate_deseq2_pca_header_ch = Channel.fromPath("$baseDir/assets/multiqc/replicate_deseq2_pca_header.txt")
replicate_deseq2_clustering_header_ch = Channel.fromPath("$baseDir/assets/multiqc/replicate_deseq2_clustering_header.txt")

sample_peak_count_header_ch = Channel.fromPath("$baseDir/assets/multiqc/sample_peak_count_header.txt")
sample_frip_score_header_ch = Channel.fromPath("$baseDir/assets/multiqc/sample_frip_score_header.txt")
sample_peak_annotation_header_ch = Channel.fromPath("$baseDir/assets/multiqc/sample_peak_annotation_header.txt")
sample_deseq2_pca_header_ch = Channel.fromPath("$baseDir/assets/multiqc/sample_deseq2_pca_header.txt")
sample_deseq2_clustering_header_ch = Channel.fromPath("$baseDir/assets/multiqc/sample_deseq2_clustering_header.txt")

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Validate inputs
if( params.fasta ){
    Channel
        .fromPath(params.fasta, checkIfExists: true)
        .ifEmpty { exit 1, "Fasta file not found: ${params.fasta}" }
        .into { fasta_bwa_index;
                fasta_genome_filter;
                fasta_markdup_metrics;
                fasta_replicate_macs_annotate;
                fasta_replicate_macs_consensus_annotate;
                fasta_sample_macs_annotate;
                fasta_sample_macs_consensus_annotate }
} else {
    exit 1, "Fasta file not specified!"
}

if( params.gtf ){
    Channel
        .fromPath(params.gtf, checkIfExists: true)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .into { gtf_makeBED12;
                gtf_replicate_macs_annotate;
                gtf_replicate_macs_consensus_annotate;
                gtf_sample_macs_annotate;
                gtf_sample_macs_consensus_annotate }
} else {
    exit 1, "GTF annotation file not specified!"
}

if( params.bwa_index ){
    Channel
        .fromPath("${params.bwa_index}/*.{amb,ann,bwt,pac,sa}", checkIfExists: true)
        .ifEmpty { exit 1, "BWA index not found: ${params.bwa_index}" }
        .into { bwa_index_read1;
                bwa_index_read2;
                bwa_index_sai_to_sam }
}

if( params.bed12 ){
    Channel
        .fromPath(params.bed12, checkIfExists: true)
        .ifEmpty { exit 1, "BED12 annotation file not found: ${params.bed12}" }
        .into { bed12_replicate_deeptools;
                bed12_sample_deeptools }
}

if ( params.blacklist ) {
    blacklist = Channel
        .fromPath(params.blacklist, checkIfExists: true)
        .ifEmpty { exit 1, "Blacklist file not found: ${params.blacklist}" }
}

/*
 * Create a channel for input read files
 */
if( params.design ){
    if ( params.singleEnd ) {
        Channel
            .fromPath(params.design, checkIfExists: true)
            .ifEmpty { exit 1, "Samples design file not found: ${params.design}" }
            .splitCsv(header:true, sep:',')
            .map { row -> [ [row.sample,"R"+row.replicate].join("_"),
                            [file(row.fastq_1)] ] }
            .groupTuple(by: [0])
            .map { group -> (1..group[1].size()).collect { value -> [ [ group[0],value ].join("_L"), group[1][value-1] ] } }
            .flatten()
            .collate( 2 )
            .map { it -> [ it[0], [ it[1] ] ] }
            .into { design_replicates_exist;
                    design_multiple_samples;
                    raw_reads_fastqc;
                    raw_reads_trimgalore }
    } else {
      Channel
          .fromPath(params.design, checkIfExists: true)
          .ifEmpty { exit 1, "Samples design file not found: ${params.design}" }
          .splitCsv(header:true, sep:',')
          .map { row -> [ [row.sample,"R"+row.replicate].join("_"),
                          [file(row.fastq_1), file(row.fastq_2)] ] }
          .groupTuple(by: [0])
          .map { group -> (1..group[1].size()).collect { value -> [ [ group[0],value ].join("_L"), group[1][value-1] ] } }
          .flatten()
          .collate( 3 )
          .map { it -> [ it[0], [ it[1], it[2] ] ] }
          .into { design_replicates_exist;
                  design_multiple_samples;
                  raw_reads_fastqc;
                  raw_reads_trimgalore }
    }
} else {
    exit 1, "Samples design file not specified!"
}

// Boolean value for replicates existing in design
replicates_exist = design_replicates_exist.map { it -> it[0][-4].toInteger() }
                                          .flatten()
                                          .max()
                                          .val > 1

// Boolean value for multiple samples existing in design
multiple_samples = design_multiple_samples.map { it -> it[0][0..-7] }
                                          .flatten()
                                          .unique()
                                          .count()
                                          .val > 1

////////////////////////////////////////////////////
/* --                   AWS                    -- */
////////////////////////////////////////////////////

// AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}

// Check workDir/outdir paths to be S3 buckets if running on AWSBatch
// related: https://github.com/nextflow-io/nextflow/issues/813
if( workflow.profile == 'awsbatch') {
    if(!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       HEADER LOG INFO                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

 nf-core/atacseq v${workflow.manifest.version}
======================================================="""

def summary = [:]
summary['Pipeline Name']          = 'nf-core/atacseq'
summary['Pipeline Version']       = workflow.manifest.version
summary['Run Name'] = custom_runName ?: workflow.runName
summary['Genome']                 = params.genome ? params.genome : 'Not supplied'
summary['Data Type']              = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Design File']            = params.design
if(params.bwa_index)  summary['BWA Index'] = params.bwa_index ? params.bwa_index : 'Not supplied'
summary['Fasta Ref']              = params.fasta
summary['GTF File']               = params.gtf
summary['BED12 File']             = params.bed12 ? params.bed12 : 'Not supplied'
if(params.blacklist) summary['Blacklist BED'] = params.blacklist
summary['Mitochondrial Contig']   = params.mito_name ? params.mito_name : 'Not supplied'
summary['MACS Genome Size']       = params.macs_gsize ? params.macs_gsize : 'Not supplied'
if(params.macs_gsize)  summary['MACS narrow peaks'] = params.narrowPeak ? 'Yes' : 'No'
if( params.skipTrimming ){
    summary['Trimming Step'] = 'Skipped'
} else {
    summary['Trim R1']    = "$params.clip_r1 bp"
    summary['Trim R2']    = "$params.clip_r2 bp"
    summary["Trim 3' R1"] = "$params.three_prime_clip_r1 bp"
    summary["Trim 3' R2"] = "$params.three_prime_clip_r2 bp"
}
summary['Fragment Size']          = "$params.fragment_size bp"
summary['Keep Mitochondrial']     = params.keepMito ? 'Yes' : 'No'
summary['Keep Duplicates']        = params.keepDups ? 'Yes' : 'No'
summary['Keep Multi-mapped']      = params.keepMultiMap ? 'Yes' : 'No'
summary['Sample-level Analysis']  = params.skipMergeBySample ? 'No' : 'Yes'
summary['Save Reference']         = params.saveReference ? 'Yes' : 'No'
summary['Save Trimmed']           = params.saveTrimmed ? 'Yes' : 'No'
summary['Save Intermeds']         = params.saveAlignedIntermediates ? 'Yes' : 'No'
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

// Show a big warning message if we're not running MACS
if (!params.macs_gsize){
    def warnstring = params.genome ? "supported for '${params.genome}'" : 'supplied'
    log.warn "=================================================================\n" +
             "  WARNING! MACS genome size parameter not $warnstring.\n" +
             "  Peak calling, annotation and differential analysis will be skipped.\n" +
             "  Please specify value for '--macs_gsize' to run these steps.\n" +
             "======================================================================="
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                     PREPARE ANNOTATION FILES                        -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * PREPROCESSING - Build BWA index
 */
if(!params.bwa_index){
    process makeBWAindex {
        tag "$fasta"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome/BWAIndex" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta_bwa_index

        output:
        file "*" into bwa_index_read1,
                      bwa_index_read2,
                      bwa_index_sai_to_sam

        script:
        """
        bwa index -a bwtsw $fasta
        """
    }
}

/*
 * PREPROCESSING - Build BED12 file
 */
if(!params.bed12){
    process makeBED12 {
        tag "$gtf"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gtf from gtf_makeBED12

        output:
        file "*.bed" into bed12_replicate_deeptools,
                          bed12_sample_deeptools

        script: // This script is bundled with the pipeline, in nf-core/atacseq/bin/
        """
        gtf2bed $gtf > ${gtf.baseName}.bed
        """
    }
}


/*
 * PREPROCESSING - Prepare genome intervals for filtering
 */
process makeGenomeFilter {
    tag "$fasta"
    publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
             saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    file fasta from fasta_genome_filter

    output:
    file "*.fai" into genome_fai                          // FAI INDEX FOR REFERENCE GENOME
    file "*.bed" into genome_filter_regions               // BED FILE WITHOUT BLACKLIST REGIONS & MITOCHONDRIAL CONTIG FOR FILTERING
    file "*.sizes" into genome_replicate_bigwig,          // CHROMOSOME SIZES FILE FOR BEDTOOLS
                        genome_sample_bigwig

    script:
    blacklist_filter = params.blacklist ? "sortBed -i ${params.blacklist} -g ${fasta}.sizes | complementBed -i stdin -g ${fasta}.sizes" : "awk '{print \$1, '0' , \$2}' OFS='\t' ${fasta}.sizes"
    name_filter = params.mito_name ? "| awk '\$1 !~ /${params.mito_name}/ {print \$0}'": ""
    mito_filter = params.keepMito ? "" : name_filter
    """
    samtools faidx $fasta
    cut -f 1,2 ${fasta}.fai > ${fasta}.sizes
    $blacklist_filter $mito_filter > ${fasta}.includable.bed
    """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                        FASTQ QC                                     -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.endsWith(".zip") ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from raw_reads_fastqc

    output:
    file "*.{zip,html}" into fastqc_reports_mqc

    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    if (params.singleEnd) {
        """
        [ ! -f  ${name}.fastq.gz ] && ln -s $reads ${name}.fastq.gz
        fastqc -q ${name}.fastq.gz
        """
    } else {
        """
        [ ! -f  ${name}_1.fastq.gz ] && ln -s ${reads[0]} ${name}_1.fastq.gz
        [ ! -f  ${name}_2.fastq.gz ] && ln -s ${reads[1]} ${name}_2.fastq.gz
        fastqc -q ${name}_1.fastq.gz
        fastqc -q ${name}_2.fastq.gz
        """
    }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                        ADAPTER TRIMMING                             -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * STEP 2 - Trim Galore!
 */
if(params.skipTrimming){
    raw_reads_trimgalore.into { trimmed_reads_aln_1;
                                trimmed_reads_aln_2;
                                trimmed_reads_sai_to_sam }
    trimgalore_results_mqc = []
    trimgalore_fastqc_reports_mqc = []
} else {
    process trim_galore {
        tag "$name"
        publishDir "${params.outdir}/trim_galore", mode: 'copy',
            saveAs: {filename ->
                if (filename.endsWith(".html")) "fastqc/$filename"
                else if (filename.endsWith(".zip")) "fastqc/zip/$filename"
                else if (filename.endsWith("trimming_report.txt")) "logs/$filename"
                else params.saveTrimmed ? filename : null
            }

        input:
        set val(name), file(reads) from raw_reads_trimgalore

        output:
        set val(name), file("*.fq.gz") into trimmed_reads_aln_1,
                                            trimmed_reads_aln_2,
                                            trimmed_reads_sai_to_sam
        file "*.txt" into trimgalore_results_mqc
        file "*.{zip,html}" into trimgalore_fastqc_reports_mqc

        script:
        // Added soft-links to original fastqs for consistent naming in MultiQC
        c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
        c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
        tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
        tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
        if (params.singleEnd) {
            """
            [ ! -f  ${name}.fastq.gz ] && ln -s $reads ${name}.fastq.gz
            trim_galore --fastqc --gzip $c_r1 $tpc_r1 ${name}.fastq.gz
            """
        } else {
            """
            [ ! -f  ${name}_1.fastq.gz ] && ln -s ${reads[0]} ${name}_1.fastq.gz
            [ ! -f  ${name}_2.fastq.gz ] && ln -s ${reads[1]} ${name}_2.fastq.gz
            trim_galore --paired --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 ${name}_1.fastq.gz ${name}_2.fastq.gz
            """
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                        ALIGN                                        -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * STEP 3.1 - align read 1 with bwa
 */
process bwa_aln_read1 {
    tag "$name"

    input:
    set val(name), file(reads) from trimmed_reads_aln_1
    file index from bwa_index_read1.collect()

    output:
    set val(name), file("*.sai") into sai_read1

    script:
    sainame = params.singleEnd ? "${name}.sai" : "${name}_1.sai"
    """
    bwa aln -t $task.cpus ${index[0].toString()[0..-5]} ${reads[0]} > $sainame
    """
}

/*
 * STEP 3.2 - align read 2 with bwa
 */
if(params.singleEnd){
    sai_to_sam = trimmed_reads_sai_to_sam.join(sai_read1, by: [0])
} else {
    process bwa_aln_read2 {
        tag "$name"

        input:
        set val(name), file(reads) from trimmed_reads_aln_2
        file index from bwa_index_read2.collect()

        output:
        set val(name), file("*.sai") into sai_read2

        script:
        """
        bwa aln -t $task.cpus ${index[0].toString()[0..-5]} ${reads[1]} > ${name}_2.sai
        """
    }
    sai_to_sam = trimmed_reads_sai_to_sam.join(sai_read1, by: [0])
                                         .join(sai_read2, by: [0])
                                         .map { it -> [it[0], it[1], [ it[2], it[3] ] ] }
}

/*
 * STEP 3.3 - convert .sai to .sam
 */
process bwa_sai_to_sam {
    tag "$name"

    input:
    set val(name), file(fastqs), file(sais) from sai_to_sam
    file index from bwa_index_sai_to_sam.collect()

    output:
    set val(name), file("*.sam") into bwa_sam

    script:
    command = params.singleEnd ? "bwa samse" : "bwa sampe"
    rg="\'@RG\\tID:${name}\\tSM:${name.toString().subSequence(0, name.length() - 3)}\\tPL:illumina\\tLB:1\\tPU:1\'"
    """
    $command -r $rg ${index[0].toString()[0..-5]} $sais $fastqs > ${name}.sam
    """
}

/*
 * STEP 3.4 - convert .sam to coordinate sorted .bam
 */
process bwa_bam {
    tag "$name"
    publishDir path: "${params.outdir}/bwa/library", mode: 'copy',
        saveAs: { filename ->
            if (params.saveAlignedIntermediates) {
                if (filename.endsWith(".flagstat")) "flagstat/$filename"
                else filename }
            }

    input:
    set val(name), file(sam) from bwa_sam

    output:
    set val(name), file("*.sorted.{bam,bam.bai}") into bwa_bam
    file "*.flagstat" into bwa_bam_flagstat

    script:
    """
    samtools view -b -h -O BAM -@ $task.cpus -o ${name}.bam $sam
    samtools sort -@ $task.cpus -o ${name}.sorted.bam -T $name ${name}.bam
    samtools index ${name}.sorted.bam
    samtools flagstat ${name}.sorted.bam > ${name}.sorted.bam.flagstat
    """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                    BAM POST-PROCESSING                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * STEP 4.1 - picard mark duplicates at library-level
 */
process markdup {
    tag "$name"
    publishDir path: "${params.outdir}/bwa/library", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith(".flagstat")) "flagstat/$filename"
            else if (filename.endsWith(".idxstats")) "idxstats/$filename"
            else if (filename.endsWith(".metrics.txt")) "picard_metrics/$filename"
            else params.saveAlignedIntermediates ? filename : null
        }

    input:
    set val(name), file(bam) from bwa_bam

    output:
    set val(name), file("*.{bam,bam.bai}") into markdup_bam_filter,
                                                markdup_bam_metrics
    file "*.{flagstat,idxstats}" into markdup_bam_stats_mqc
    file "*.txt" into markdup_metrics_mqc

    script:
    prefix="${name}.mkD"
    if( !task.memory ){
        log.info "[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this."
        avail_mem = 3
    } else {
        avail_mem = task.memory.toGiga()
    }
    """
    picard -Xmx${avail_mem}g MarkDuplicates \\
           INPUT=${bam[0]} \\
           OUTPUT=${prefix}.sorted.bam \\
           ASSUME_SORTED=true \\
           REMOVE_DUPLICATES=false \\
           METRICS_FILE=${prefix}.MarkDuplicates.metrics.txt \\
           VALIDATION_STRINGENCY=LENIENT \\
           TMP_DIR=tmp

    samtools index ${prefix}.sorted.bam
    samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
    samtools idxstats ${prefix}.sorted.bam > ${prefix}.sorted.bam.idxstats
    """
}

/*
 * STEP 4.2 picard collectmultiplemetrics at library-level
 */
process markdup_collectmetrics {
    tag "$name"
    publishDir path: "${params.outdir}/bwa/library", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith("_metrics")) "picard_metrics/$filename"
            else if (filename.endsWith(".pdf")) "picard_metrics/pdf/$filename"
            else null
        }

    input:
    set val(name), file(bam) from markdup_bam_metrics
    file fasta from fasta_markdup_metrics.collect()

    output:
    file "*_metrics" into markdup_collectmetrics_mqc
    file "*.pdf" into markdup_collectmetrics_pdf

    script:
    prefix="${name}.mkD"
    if( !task.memory ){
        log.info "[Picard CollectMultipleMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this."
        avail_mem = 3
    } else {
        avail_mem = task.memory.toGiga()
    }
    """
    picard -Xmx${avail_mem}g CollectMultipleMetrics \\
           INPUT=${bam[0]} \\
           OUTPUT=${prefix}.CollectMultipleMetrics \\
           REFERENCE_SEQUENCE=$fasta \\
           VALIDATION_STRINGENCY=LENIENT \\
           TMP_DIR=tmp
    """
}

/*
 * STEP 4.3 filter bam
 */
process filter_bam {
    tag "$name"
    publishDir path: "${params.outdir}/bwa/library", mode: 'copy',
        saveAs: { filename ->
            if (params.singleEnd || params.saveAlignedIntermediates) {
                if (filename.endsWith(".flagstat")) "flagstat/$filename"
                else filename }
            }

    input:
    set val(name), file(bam) from markdup_bam_filter
    file bed from genome_filter_regions.collect()
    file bamtools_filter_se_config from bamtools_filter_se_config_ch.collect()
    file bamtools_filter_pe_config from bamtools_filter_pe_config_ch.collect()

    output:
    set val(name), file("*.{bam,bam.bai}") into filter_bam
    file "*.flagstat" into filter_bam_flagstat

    script:
    prefix = params.singleEnd ? "${name}.clN" : "${name}.flT"
    filter_params = params.singleEnd ? "-F 0x004 -F 0x0100" : "-f 0x001 -F 0x004 -F 0x0008 -F 0x0100"
    dup_params = params.keepDups ? "" : "-F 0x0400"
    multimap_params = params.keepMultiMap ? "" : "-q 1"
    bamtools_filter_config = params.singleEnd ? bamtools_filter_se_config : bamtools_filter_pe_config
    name_sort_bam = params.singleEnd ? "" : "samtools sort -n -@ $task.cpus -o ${prefix}.bam -T $prefix ${prefix}.sorted.bam"
    """
    samtools view \\
             $filter_params \\
             $dup_params \\
             $multimap_params \\
             -L $bed \\
             -b ${bam[0]} \\
             | bamtools filter \\
                        -out ${prefix}.sorted.bam \\
                        -script $bamtools_filter_config

    samtools index ${prefix}.sorted.bam
    samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
    $name_sort_bam
    """
}


/*
 * STEP 4.4 remove orphan reads from paired-end BAM
 */
if(params.singleEnd){
    filter_bam.into { rm_orphan_bam_replicate;
                      rm_orphan_bam_sample }
    rm_orphan_flagstat_mqc = filter_bam_flagstat
} else {
    process rm_orphan {
        tag "$name"
        publishDir path: "${params.outdir}/bwa/library", mode: 'copy',
            saveAs: { filename ->
                if (filename.endsWith(".flagstat")) "flagstat/$filename"
                else filename
            }

        input:
        set val(name), file(bam) from filter_bam

        output:
        set val(name), file("*.sorted.{bam,bam.bai}") into rm_orphan_bam_replicate,
                                                           rm_orphan_bam_sample
        file "*.flagstat" into rm_orphan_flagstat_mqc

        script: // This script is bundled with the pipeline, in nf-core/atacseq/bin/
        prefix="${name}.clN"
        """
        bampe_rm_orphan.py ${bam[0]} ${prefix}.bam --only_fr_pairs
        samtools sort -@ $task.cpus -o ${prefix}.sorted.bam -T $prefix ${prefix}.bam
        samtools index ${prefix}.sorted.bam
        samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
        """
    }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                    MERGE REPLICATE BAM                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * STEP 5.1 merge bam files for all libraries from same replicate
 */
rm_orphan_bam_replicate.map { it -> [ it[0].toString().subSequence(0, it[0].length() - 3), it[1] ] }
                       .groupTuple(by: [0])
                       .map { it ->  [ it[0], it[1].flatten() ] }
                       .set { rm_orphan_bam_replicate }

process merge_replicate {
    tag "$name"
    publishDir "${params.outdir}/bwa/replicate", mode: 'copy',
        saveAs: {filename ->
                    if (filename.endsWith(".flagstat")) "flagstat/$filename"
                    else if (filename.endsWith(".metrics.txt")) "picard_metrics/$filename"
                    else filename
                }

    input:
    set val(name), file(bams) from rm_orphan_bam_replicate

    output:
    set val(name), file("*${prefix}.sorted.{bam,bam.bai}") into merge_replicate_bam,
                                                                merge_replicate_bam_bigwig,
                                                                merge_replicate_bam_macs
    set val(name), file("*.flagstat") into merge_replicate_flagstat_mqc,
                                           merge_replicate_flagstat_bigwig,
                                           merge_replicate_flagstat_macs
    file "*.txt" into merge_replicate_metrics_mqc

    script:
    prefix="${name}.mRp"
    bam_files = bams.findAll { it.toString().endsWith('.bam') }.sort()
    rmdup = params.keepDups ? "false" : "true"
    if( !task.memory ){
        log.info "[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this."
        avail_mem = 3
    } else {
        avail_mem = task.memory.toGiga()
    }
    if (bam_files.size() > 1) {
        """
        picard -Xmx${avail_mem}g MergeSamFiles \\
               ${'INPUT='+bam_files.join(' INPUT=')} \\
               OUTPUT=${name}.sorted.bam \\
               SORT_ORDER=coordinate \\
               VALIDATION_STRINGENCY=LENIENT \\
               TMP_DIR=tmp
        samtools index ${name}.sorted.bam

        picard -Xmx${avail_mem}g MarkDuplicates \\
               INPUT=${name}.sorted.bam \\
               OUTPUT=${prefix}.sorted.bam \\
               ASSUME_SORTED=true \\
               REMOVE_DUPLICATES=$rmdup \\
               METRICS_FILE=${prefix}.MarkDuplicates.metrics.txt \\
               VALIDATION_STRINGENCY=LENIENT \\
               TMP_DIR=tmp

        samtools index ${prefix}.sorted.bam
        samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
        """
    } else {
      """
      ln -s ${bams[0]} ${prefix}.sorted.bam
      ln -s ${bams[1]} ${prefix}.sorted.bam.bai
      touch ${prefix}.MarkDuplicates.metrics.txt
      samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
      """
    }
}

/*
 * STEP 5.2 sort paired-end merged bam file by name for featurecounts
 */
if(params.singleEnd){
    merge_replicate_bam.into { replicate_name_bam_replicate_counts;
                               replicate_name_bam_sample_counts }
} else {
    process replicate_name_bam {
        tag "$name"

        input:
        set val(name), file(bam) from merge_replicate_bam

        output:
        set val(name), file("*.bam") into replicate_name_bam_replicate_counts,
                                          replicate_name_bam_sample_counts

        script:
        prefix="${name}.mRp"
        """
        samtools sort -n -@ $task.cpus -o ${prefix}.bam -T $prefix ${bam[0]}
        """
    }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                 MERGE REPLICATE BAM POST-ANALYSIS                   -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * STEP 5.3 Read depth normalised bigWig
 */
process replicate_bigwig {
    tag "$name"
    publishDir "${params.outdir}/bwa/replicate/bigwig", mode: 'copy',
        saveAs: {filename ->
                    if (filename.endsWith(".txt")) "scale/$filename"
                    else if (filename.endsWith(".bigWig")) "$filename"
                    else null
                }

    input:
    set val(name), file(bam), file(flagstat) from merge_replicate_bam_bigwig.join(merge_replicate_flagstat_bigwig, by: [0])
    file sizes from genome_replicate_bigwig.collect()

    output:
    file "*.bigWig" into replicate_bigwig,
                         replicate_bigwig_igv
    file "*.txt" into replicate_bigwig_scale

    script:
    prefix="${name}.mRp"
    pe_fragment = params.singleEnd ? "" : "-pc"
    extend = (params.singleEnd && params.fragment_size > 0) ? "-fs ${params.fragment_size}" : ''
    """
    SCALE_FACTOR=\$(grep 'mapped (' $flagstat | awk '{print 1000000/\$1}')
    echo \$SCALE_FACTOR > ${prefix}.scale_factor.txt
    genomeCoverageBed -ibam ${bam[0]} -bg -trackline -scale \$SCALE_FACTOR $pe_fragment $extend >  ${prefix}.bedGraph
    wigToBigWig -clip ${prefix}.bedGraph $sizes ${prefix}.bigWig
    """
}

/*
 * STEP 5.4 generate TSS profiles with deeptools
 */
process replicate_tss_plot {
    publishDir "${params.outdir}/bwa/replicate/deeptools", mode: 'copy'

    input:
    file bigwigs from replicate_bigwig.collect()
    file bed from bed12_replicate_deeptools.collect()

    output:
    file "*.png" into replicate_tss_plot

    script:
    suffix='mRp'
    prefix="plotProfile.${suffix}"
    """
    computeMatrix reference-point \\
                  -R $bed \\
                  -S ${bigwigs.collect{it.toString()}.sort().join(' ')} \\
                  --samplesLabel ${bigwigs.collect{it.toString()}.sort().join(' ').replaceAll(".${suffix}.bigWig","")} \\
                  -out ${prefix}.TSS.mat.gz \\
                  --referencePoint TSS \\
                  -a 3000 \\
                  -b 3000 \\
                  -p $task.cpus

    plotProfile -m ${prefix}.TSS.mat.gz \\
                -out ${prefix}.TSS.png
    """
}

/*
 * STEP 5.5.1 Call peaks with MACS2 and calculate FRiP score
 */
process replicate_macs {
    tag "$name"
    publishDir "${params.outdir}/bwa/replicate/macs", mode: 'copy',
        saveAs: {filename ->
                    if (filename.endsWith(".tsv")) "qc/$filename"
                    else filename
                }

    input:
    set val(name), file(bam), file(flagstat) from merge_replicate_bam_macs.join(merge_replicate_flagstat_macs, by: [0])
    file replicate_peak_count_header from replicate_peak_count_header_ch.collect()
    file replicate_frip_score_header from replicate_frip_score_header_ch.collect()

    output:
    file "*.{bed,xls,gappedPeak}" into replicate_macs_output
    set val(name), file("*$peakext") into replicate_macs_peaks_homer,
                                          replicate_macs_peaks_qc,
                                          replicate_macs_peaks_consensus,
                                          replicate_macs_peaks_igv
    file "*_mqc.tsv" into replicate_macs_peak_mqc

    when: params.macs_gsize

    script:
    prefix="${name}.mRp"
    peakext = params.narrowPeak ? ".narrowPeak" : ".broadPeak"
    broad = params.narrowPeak ? '' : "--broad"
    format = params.singleEnd ? "BAM" : "BAMPE"
    """
    macs2 callpeak \\
         -t ${bam[0]} \\
         $broad \\
         -f $format \\
         -g ${params.macs_gsize} \\
         -n ${prefix} \\
         --keep-dup all \\
         --nomodel

    cat ${prefix}_peaks${peakext} | wc -l | awk -v OFS='\t' '{ print "${name}", \$1 }' | cat $replicate_peak_count_header - > ${prefix}_peaks.count_mqc.tsv

    READS_IN_PEAKS=\$(intersectBed -a ${bam[0]} -b ${prefix}_peaks${peakext} -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
    grep 'mapped (' $flagstat | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${name}", a/\$1}' | cat $replicate_frip_score_header - > ${prefix}_peaks.FRiP_mqc.tsv
    """
}

/*
 * STEP 5.5.2 annotate peaks with homer
 */
process replicate_macs_annotate {
    tag "$name"
    publishDir "${params.outdir}/bwa/replicate/macs", mode: 'copy'

    input:
    set val(name), file(peak) from replicate_macs_peaks_homer
    file fasta from fasta_replicate_macs_annotate.collect()
    file gtf from gtf_replicate_macs_annotate.collect()

    output:
    file "*.txt" into replicate_macs_annotate

    when: params.macs_gsize

    script:
    prefix="${name}.mRp"
    """
    annotatePeaks.pl $peak \\
                     $fasta \\
                     -gid \\
                     -gtf $gtf \\
                     > ${prefix}_peaks.annotatePeaks.txt
    """
}

/*
 * STEP 5.5.3 aggregated qc plots for peaks, frip and annotation
 */
process replicate_macs_qc {
   publishDir "${params.outdir}/bwa/replicate/macs/qc", mode: 'copy'

   input:
   file peaks from replicate_macs_peaks_qc.collect{ it[1] }
   file annos from replicate_macs_annotate.collect()
   file replicate_peak_annotation_header from replicate_peak_annotation_header_ch.collect()

   output:
   file "*.{txt,pdf}" into replicate_macs_qc
   file "*.tsv" into replicate_macs_qc_mqc

   when: params.macs_gsize

   script:  // This script is bundled with the pipeline, in nf-core/atacseq/bin/
   suffix='mRp'
   peakext = params.narrowPeak ? ".narrowPeak" : ".broadPeak"
   """
   plot_macs_qc.r -i ${peaks.join(',')} \\
                  -s ${peaks.join(',').replaceAll(".${suffix}_peaks${peakext}","")} \\
                  -o ./ \\
                  -p macs_peak.${suffix}

   plot_homer_annotatepeaks.r -i ${annos.join(',')} \\
                              -s ${annos.join(',').replaceAll(".${suffix}_peaks.annotatePeaks.txt","")} \\
                              -o ./ \\
                              -p macs_annotatePeaks.${suffix}

   cat $replicate_peak_annotation_header macs_annotatePeaks.${suffix}.summary.txt > macs_annotatePeaks.${suffix}.summary_mqc.tsv
   """
}

/*
 * STEP 5.5.4 consensus peaks across samples, create boolean filtering file, saf file for featurecounts and UpSetR plot for intersection
 */
process replicate_macs_consensus {
    publishDir "${params.outdir}/bwa/replicate/macs/consensus", mode: 'copy'

    input:
    file peaks from replicate_macs_peaks_consensus.collect{ it[1] }

    output:
    file "*.bed" into replicate_macs_consensus_bed,
                      replicate_macs_consensus_igv
    file "*.boolean.txt" into replicate_macs_consensus_bool
    file "*.saf" into replicate_macs_consensus_saf
    file "*.intersect.{txt,plot.pdf}" into replicate_macs_consensus_intersect

    when: params.macs_gsize && (multiple_samples || replicates_exist)

    script: // scripts are bundled with the pipeline, in nf-core/atacseq/bin/
    suffix='mRp'
    prefix="consensus_peaks.${suffix}"
    peakext = params.narrowPeak ? ".narrowPeak" : ".broadPeak"
    mergecols = params.narrowPeak ? (2..10).join(',') : (2..9).join(',')
    collapsecols = params.narrowPeak ? (["collapse"]*9).join(',') : (["collapse"]*8).join(',')
    expandparam = params.narrowPeak ? "--is_narrow_peak" : ""
    """
    sort -k1,1 -k2,2n ${peaks.collect{it.toString()}.sort().join(' ')} \\
         | mergeBed -c $mergecols -o $collapsecols > ${prefix}.txt

    macs2_merged_expand.py ${prefix}.txt \\
                           ${peaks.collect{it.toString()}.sort().join(',').replaceAll("_peaks${peakext}","")} \\
                           ${prefix}.boolean.txt \\
                           --min_samples 1 \\
                           $expandparam

    awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$2, \$3, \$4, "0", "+" }' ${prefix}.boolean.txt > ${prefix}.bed

    echo -e "GeneID\tChr\tStart\tEnd\tStrand" > ${prefix}.saf
    awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$4, \$1, \$2, \$3,  "+" }' ${prefix}.boolean.txt >> ${prefix}.saf

    sed -i 's/.${suffix}//g' ${prefix}.boolean.intersect.txt
    plot_peak_intersect.r -i ${prefix}.boolean.intersect.txt -o ${prefix}.boolean.intersect.plot.pdf
    """
}

/*
 * STEP 5.5.5 annotate consensus peaks with homer, and add annotation to boolean output file
 */
process replicate_macs_consensus_annotate {
    publishDir "${params.outdir}/bwa/replicate/macs/consensus", mode: 'copy'

    input:
    file bed from replicate_macs_consensus_bed
    file bool from replicate_macs_consensus_bool
    file fasta from fasta_replicate_macs_consensus_annotate
    file gtf from gtf_replicate_macs_consensus_annotate

    output:
    file "*.annotatePeaks.txt" into replicate_macs_consensus_annotate

    when: params.macs_gsize && (multiple_samples || replicates_exist)

    script:
    prefix="consensus_peaks.mRp"
    """
    annotatePeaks.pl $bed \\
                     $fasta \\
                     -gid \\
                     -gtf $gtf \\
                     > ${prefix}.annotatePeaks.txt

    cut -f2- ${prefix}.annotatePeaks.txt | awk 'NR==1; NR > 1 {print \$0 | "sort -k1,1 -k2,2n"}' | cut -f6- > tmp.txt
    paste $bool tmp.txt > ${prefix}.boolean.annotatePeaks.txt
    """
}

/*
 * STEP 5.5.6 count reads in consensus peaks with featurecounts
 */
process replicate_macs_consensus_deseq {
    publishDir "${params.outdir}/bwa/replicate/macs/consensus/deseq2", mode: 'copy'

    input:
    file bams from replicate_name_bam_replicate_counts.collect{ it[1] }
    file saf from replicate_macs_consensus_saf.collect()
    file replicate_deseq2_pca_header from replicate_deseq2_pca_header_ch.collect()
    file replicate_deseq2_clustering_header from replicate_deseq2_clustering_header_ch.collect()

    output:
    file "*featureCounts.txt" into replicate_macs_consensus_counts
    file "*featureCounts.txt.summary" into replicate_macs_consensus_counts_mqc
    file "*.{RData,results.txt,pdf,log}" into replicate_macs_consensus_deseq_results
    file "sizeFactors" into replicate_macs_consensus_deseq_factors
    file "*vs*/*.{pdf,txt}" into replicate_macs_consensus_deseq_comp_results
    file "*vs*/*.bed" into replicate_macs_consensus_deseq_comp_bed
    file "*.tsv" into replicate_macs_consensus_deseq_mqc

    when: params.macs_gsize && multiple_samples && replicates_exist

    script:
    prefix="consensus_peaks.mRp"
    bam_files = bams.findAll { it.toString().endsWith('.bam') }.sort()
    bam_ext = params.singleEnd ? ".mRp.sorted.bam" : ".mRp.bam"
    pe_params = params.singleEnd ? '' : "-p --donotsort"
    """
    featureCounts -F SAF \\
                  -O \\
                  --fracOverlap 0.2 \\
                  -T $task.cpus \\
                  $pe_params \\
                  -a $saf \\
                  -o ${prefix}.featureCounts.txt \\
                  ${bam_files.join(' ')}

    featurecounts_deseq2.r -i ${prefix}.featureCounts.txt -b '$bam_ext' -o ./ -p $prefix -s .mRp

    cat $replicate_deseq2_pca_header ${prefix}.pca.vals.txt > ${prefix}.pca.vals_mqc.tsv
    cat $replicate_deseq2_clustering_header ${prefix}.sample.dists.txt > ${prefix}.sample.dists_mqc.tsv
    """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                    MERGE SAMPLE BAM                                 -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * STEP 6.1 merge bam files for all libraries from same sample
 */
rm_orphan_bam_sample.map { it -> [ it[0].toString().subSequence(0, it[0].length() - 6), it[1] ] }
                    .groupTuple(by: [0])
                    .map { it ->  [ it[0], it[1].flatten() ] }
                    .set { rm_orphan_bam_sample }

process merge_sample {
    tag "$name"
    publishDir "${params.outdir}/bwa/sample", mode: 'copy',
        saveAs: {filename ->
                    if (filename.endsWith(".flagstat")) "flagstat/$filename"
                    else if (filename.endsWith(".metrics.txt")) "picard_metrics/$filename"
                    else filename
                }

    input:
    set val(name), file(bams) from rm_orphan_bam_sample

    output:
    set val(name), file("*${prefix}.sorted.{bam,bam.bai}") into merge_sample_bam_bigwig,
                                                                merge_sample_bam_macs
    set val(name), file("*.flagstat") into merge_sample_flagstat_mqc,
                                           merge_sample_flagstat_bigwig,
                                           merge_sample_flagstat_macs
    file "*.txt" into merge_sample_metrics_mqc

    when: !skipMergeBySample && replicates_exist

    script:
    prefix="${name}.mSm"
    bam_files = bams.findAll { it.toString().endsWith('.bam') }.sort()
    rmdup = params.keepDups ? "false" : "true"
    if( !task.memory ){
        log.info "[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this."
        avail_mem = 3
    } else {
        avail_mem = task.memory.toGiga()
    }
    if (bam_files.size() > 1) {
        """
        picard -Xmx${avail_mem}g MergeSamFiles \\
               ${'INPUT='+bam_files.join(' INPUT=')} \\
               OUTPUT=${name}.sorted.bam \\
               SORT_ORDER=coordinate \\
               VALIDATION_STRINGENCY=LENIENT \\
               TMP_DIR=tmp
        samtools index ${name}.sorted.bam

        picard -Xmx${avail_mem}g MarkDuplicates \\
               INPUT=${name}.sorted.bam \\
               OUTPUT=${prefix}.sorted.bam \\
               ASSUME_SORTED=true \\
               REMOVE_DUPLICATES=$rmdup \\
               METRICS_FILE=${prefix}.MarkDuplicates.metrics.txt \\
               VALIDATION_STRINGENCY=LENIENT \\
               TMP_DIR=tmp

        samtools index ${prefix}.sorted.bam
        samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
        """
    } else {
      """
      ln -s ${bams[0]} ${prefix}.sorted.bam
      ln -s ${bams[1]} ${prefix}.sorted.bam.bai
      touch ${prefix}.MarkDuplicates.metrics.txt
      samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
      """
    }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                 MERGE SAMPLE BAM POST-ANALYSIS                      -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * STEP 6.2 Read depth normalised bigWig
 */
process sample_bigwig {
    tag "$name"
    publishDir "${params.outdir}/bwa/sample/bigwig", mode: 'copy',
        saveAs: {filename ->
                    if (filename.endsWith(".txt")) "scale/$filename"
                    else if (filename.endsWith(".bigWig")) "$filename"
                    else null
                }

    input:
    set val(name), file(bam), file(flagstat) from merge_sample_bam_bigwig.join(merge_sample_flagstat_bigwig, by: [0])
    file sizes from genome_sample_bigwig.collect()

    output:
    file "*.bigWig" into sample_bigwig,
                         sample_bigwig_igv
    file "*.txt" into sample_bigwig_scale

    when: !skipMergeBySample && replicates_exist

    script:
    prefix="${name}.mSm"
    pe_fragment = params.singleEnd ? "" : "-pc"
    extend = (params.singleEnd && params.fragment_size > 0) ? "-fs ${params.fragment_size}" : ''
    """
    SCALE_FACTOR=\$(grep 'mapped (' $flagstat | awk '{print 1000000/\$1}')
    echo \$SCALE_FACTOR > ${prefix}.scale_factor.txt
    genomeCoverageBed -ibam ${bam[0]} -bg -trackline -scale \$SCALE_FACTOR $pe_fragment $extend >  ${prefix}.bedGraph
    wigToBigWig -clip ${prefix}.bedGraph $sizes ${prefix}.bigWig
    """
}

/*
 * STEP 6.3 generate TSS profiles with deeptools
 */
process sample_tss_plot {
    publishDir "${params.outdir}/bwa/sample/deeptools", mode: 'copy'

    input:
    file bigwigs from sample_bigwig.collect()
    file bed from bed12_sample_deeptools.collect()

    output:
    file "*.png" into sample_tss_plot

    when: !skipMergeBySample && replicates_exist

    script:
    suffix='mSm'
    prefix="plotProfile.${suffix}"
    """
    computeMatrix reference-point \\
                  -R $bed \\
                  -S ${bigwigs.collect{it.toString()}.sort().join(' ')} \\
                  --samplesLabel ${bigwigs.collect{it.toString()}.sort().join(' ').replaceAll(".${suffix}.bigWig","")} \\
                  -out ${prefix}.TSS.mat.gz \\
                  --referencePoint TSS \\
                  -a 3000 \\
                  -b 3000 \\
                  -p $task.cpus

    plotProfile -m ${prefix}.TSS.mat.gz \\
                -out ${prefix}.TSS.png
    """
}

/*
 * STEP 6.4.1 Call peaks with MACS2 and calculate FRiP score
 */
process sample_macs {
    tag "$name"
    publishDir "${params.outdir}/bwa/sample/macs", mode: 'copy',
        saveAs: {filename ->
                    if (filename.endsWith(".tsv")) "qc/$filename"
                    else filename
                }

    input:
    set val(name), file(bam), file(flagstat) from merge_sample_bam_macs.join(merge_sample_flagstat_macs, by: [0])
    file sample_peak_count_header from sample_peak_count_header_ch.collect()
    file sample_frip_score_header from sample_frip_score_header_ch.collect()

    output:
    file "*.{bed,xls,gappedPeak}" into sample_macs_output
    set val(name), file("*$peakext") into sample_macs_peaks_homer,
                                          sample_macs_peaks_qc,
                                          sample_macs_peaks_consensus,
                                          sample_macs_peaks_igv
    file "*_mqc.tsv" into sample_macs_peak_mqc

    when: !skipMergeBySample && replicates_exist && params.macs_gsize

    script:
    prefix="${name}.mSm"
    peakext = params.narrowPeak ? ".narrowPeak" : ".broadPeak"
    broad = params.narrowPeak ? '' : "--broad"
    format = params.singleEnd ? "BAM" : "BAMPE"
    """
    macs2 callpeak \\
         -t ${bam[0]} \\
         $broad \\
         -f $format \\
         -g ${params.macs_gsize} \\
         -n ${prefix} \\
         --keep-dup all \\
         --nomodel

    cat ${prefix}_peaks${peakext} | wc -l | awk -v OFS='\t' '{ print "${name}", \$1 }' | cat $sample_peak_count_header - > ${prefix}_peaks.count_mqc.tsv

    READS_IN_PEAKS=\$(intersectBed -a ${bam[0]} -b ${prefix}_peaks${peakext} -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
    grep 'mapped (' $flagstat | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${name}", a/\$1}' | cat $sample_frip_score_header - > ${prefix}_peaks.FRiP_mqc.tsv
    """
}

/*
 * STEP 6.4.2 annotate peaks with homer
 */
process sample_macs_annotate {
    tag "$name"
    publishDir "${params.outdir}/bwa/sample/macs", mode: 'copy'

    input:
    set val(name), file(peak) from sample_macs_peaks_homer
    file fasta from fasta_sample_macs_annotate.collect()
    file gtf from gtf_sample_macs_annotate.collect()

    output:
    file "*.txt" into sample_macs_annotate

    when: !skipMergeBySample && replicates_exist && params.macs_gsize

    script:
    prefix="${name}.mSm"
    """
    annotatePeaks.pl $peak \\
                     $fasta \\
                     -gid \\
                     -gtf $gtf \\
                     > ${prefix}_peaks.annotatePeaks.txt
    """
}

/*
 * STEP 6.4.3 aggregated qc plots for peaks, frip and annotation
 */
process sample_macs_qc {
   publishDir "${params.outdir}/bwa/sample/macs/qc", mode: 'copy'

   input:
   file peaks from sample_macs_peaks_qc.collect{ it[1] }
   file annos from sample_macs_annotate.collect()
   file sample_peak_annotation_header from sample_peak_annotation_header_ch.collect()

   output:
   file "*.{txt,pdf}" into sample_macs_qc
   file "*.tsv" into sample_macs_qc_mqc

   when: !skipMergeBySample && replicates_exist && params.macs_gsize

   script:  // This script is bundled with the pipeline, in nf-core/atacseq/bin/
   suffix='mSm'
   peakext = params.narrowPeak ? ".narrowPeak" : ".broadPeak"
   """
   plot_macs_qc.r -i ${peaks.join(',')} \\
                  -s ${peaks.join(',').replaceAll(".${suffix}_peaks${peakext}","")} \\
                  -o ./ \\
                  -p macs_peak.${suffix}

   plot_homer_annotatepeaks.r -i ${annos.join(',')} \\
                              -s ${annos.join(',').replaceAll(".${suffix}_peaks.annotatePeaks.txt","")} \\
                              -o ./ \\
                              -p macs_annotatePeaks.${suffix}

   cat $sample_peak_annotation_header macs_annotatePeaks.${suffix}.summary.txt > macs_annotatePeaks.${suffix}.summary_mqc.tsv
   """
}

/*
 * STEP 6.4.4 consensus peaks across samples, create boolean filtering file, saf file for featurecounts and UpSetR plot for intersection
 */
process sample_macs_consensus {
    publishDir "${params.outdir}/bwa/sample/macs/consensus", mode: 'copy'

    input:
    file peaks from sample_macs_peaks_consensus.collect{ it[1] }

    output:
    file "*.bed" into sample_macs_consensus_bed,
                      sample_macs_consensus_igv
    file "*.boolean.txt" into sample_macs_consensus_bool
    file "*.saf" into sample_macs_consensus_saf
    file "*.intersect.{txt,plot.pdf}" into sample_macs_consensus_intersect

    when: !skipMergeBySample && replicates_exist && params.macs_gsize && multiple_samples

    script: // scripts are bundled with the pipeline, in nf-core/atacseq/bin/
    suffix='mSm'
    prefix="consensus_peaks.${suffix}"
    peakext = params.narrowPeak ? ".narrowPeak" : ".broadPeak"
    mergecols = params.narrowPeak ? (2..10).join(',') : (2..9).join(',')
    collapsecols = params.narrowPeak ? (["collapse"]*9).join(',') : (["collapse"]*8).join(',')
    expandparam = params.narrowPeak ? "--is_narrow_peak" : ""
    """
    sort -k1,1 -k2,2n ${peaks.collect{it.toString()}.sort().join(' ')} \\
         | mergeBed -c $mergecols -o $collapsecols > ${prefix}.txt

    macs2_merged_expand.py ${prefix}.txt \\
                           ${peaks.collect{it.toString()}.sort().join(',').replaceAll("_peaks${peakext}","")} \\
                           ${prefix}.boolean.txt \\
                           --min_samples 1 \\
                           $expandparam

    awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$2, \$3, \$4, "0", "+" }' ${prefix}.boolean.txt > ${prefix}.bed

    echo -e "GeneID\tChr\tStart\tEnd\tStrand" > ${prefix}.saf
    awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$4, \$1, \$2, \$3,  "+" }' ${prefix}.boolean.txt >> ${prefix}.saf

    sed -i 's/.${suffix}//g' ${prefix}.boolean.intersect.txt
    plot_peak_intersect.r -i ${prefix}.boolean.intersect.txt -o ${prefix}.boolean.intersect.plot.pdf
    """
}

/*
 * STEP 6.4.5 annotate consensus peaks with homer, and add annotation to boolean output file
 */
process sample_macs_consensus_annotate {
    publishDir "${params.outdir}/bwa/sample/macs/consensus", mode: 'copy'

    input:
    file bed from sample_macs_consensus_bed
    file bool from sample_macs_consensus_bool
    file fasta from fasta_sample_macs_consensus_annotate
    file gtf from gtf_sample_macs_consensus_annotate

    output:
    file "*.annotatePeaks.txt" into sample_macs_consensus_annotate

    when: !skipMergeBySample && replicates_exist && params.macs_gsize && multiple_samples

    script:
    prefix="consensus_peaks.mSm"
    """
    annotatePeaks.pl $bed \\
                     $fasta \\
                     -gid \\
                     -gtf $gtf \\
                     > ${prefix}.annotatePeaks.txt

    cut -f2- ${prefix}.annotatePeaks.txt | awk 'NR==1; NR > 1 {print \$0 | "sort -k1,1 -k2,2n"}' | cut -f6- > tmp.txt
    paste $bool tmp.txt > ${prefix}.boolean.annotatePeaks.txt
    """
}

/*
 * STEP 6.4.6 count reads in consensus peaks with featurecounts
 */
process sample_macs_consensus_deseq {
    publishDir "${params.outdir}/bwa/sample/macs/consensus/deseq2", mode: 'copy'

    input:
    file bams from replicate_name_bam_sample_counts.collect{ it[1] }
    file saf from sample_macs_consensus_saf.collect()
    file sample_deseq2_pca_header from sample_deseq2_pca_header_ch.collect()
    file sample_deseq2_clustering_header from sample_deseq2_clustering_header_ch.collect()

    output:
    file "*featureCounts.txt" into sample_macs_consensus_counts
    file "*featureCounts.txt.summary" into sample_macs_consensus_counts_mqc
    file "*.{RData,results.txt,pdf,log}" into sample_macs_consensus_deseq_results
    file "sizeFactors" into sample_macs_consensus_deseq_factors
    file "*vs*/*.{pdf,txt}" into sample_macs_consensus_deseq_comp_results
    file "*vs*/*.bed" into sample_macs_consensus_deseq_comp_bed
    file "*.tsv" into sample_macs_consensus_deseq_mqc

    when: !skipMergeBySample && replicates_exist && params.macs_gsize && multiple_samples

    script:
    prefix="consensus_peaks.mSm"
    bam_files = bams.findAll { it.toString().endsWith('.bam') }.sort()
    bam_ext = params.singleEnd ? ".mRp.sorted.bam" : ".mRp.bam"
    pe_params = params.singleEnd ? '' : "-p --donotsort"
    """
    featureCounts -F SAF \\
                  -O \\
                  --fracOverlap 0.2 \\
                  -T $task.cpus \\
                  $pe_params \\
                  -a $saf \\
                  -o ${prefix}.featureCounts.txt \\
                  ${bam_files.join(' ')}

    featurecounts_deseq2.r -i ${prefix}.featureCounts.txt -b '$bam_ext' -o ./ -p $prefix -s .mSm

    cat $sample_deseq2_pca_header ${prefix}.pca.vals.txt > ${prefix}.pca.vals_mqc.tsv
    cat $sample_deseq2_clustering_header ${prefix}.sample.dists.txt > ${prefix}.sample.dists_mqc.tsv
    """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                          MULTIQC                                    -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

def create_workflow_summary(summary) {

    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-atacseq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/atacseq Workflow Summary'
    section_href: 'https://github.com/nf-core/atacseq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file "software_versions_mqc.yaml" into software_versions_yaml

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    trim_galore --version > v_trim_galore.txt
    echo \$(bwa 2>&1) > v_bwa.txt
    samtools --version > v_samtools.txt
    bedtools --version > v_bedtools.txt
    echo \$(bamtools --version 2>&1) > v_bamtools.txt
    picard MarkDuplicates --version &> v_picard.txt  || true
    echo \$(R --version 2>&1) > v_R.txt
    python -c "import pysam; print(pysam.__version__)" > v_pysam.txt
    echo \$(macs2 --version 2>&1) > v_macs2.txt
    touch v_homer.txt
    echo \$(featureCounts -v 2>&1) > v_featurecounts.txt
    echo \$(computeMatrix --version 2>&1) > v_deeptools.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

/*
 * STEP 7 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    file multiqc_config from multiqc_config_ch.collect()
    file ('fastqc/*') from fastqc_reports_mqc.collect()
    file ('trimgalore/*') from trimgalore_results_mqc.collect()
    file ('trimgalore/fastqc/*') from trimgalore_fastqc_reports_mqc.collect()
    file ('alignment/library/*') from markdup_bam_stats_mqc.collect()
    file ('alignment/library/picard_metrics/*') from markdup_metrics_mqc.collect()
    file ('alignment/library/picard_metrics/*') from markdup_collectmetrics_mqc.collect()
    file ('alignment/library/*') from rm_orphan_flagstat_mqc.collect()
    file ('alignment/replicate/*') from merge_replicate_flagstat_mqc.collect{it[1]}
    file ('alignment/replicate/picard_metrics/*') from merge_replicate_metrics_mqc.collect()
    file ('macs/replicate/*') from replicate_macs_peak_mqc.collect().ifEmpty([])
    file ('macs/replicate/*') from replicate_macs_qc_mqc.collect().ifEmpty([])
    file ('macs/replicate/consensus/*') from replicate_macs_consensus_counts_mqc.collect().ifEmpty([])
    file ('macs/replicate/consensus/*') from replicate_macs_consensus_deseq_mqc.collect().ifEmpty([])
    file ('alignment/sample/*') from merge_sample_flagstat_mqc.collect{it[1]}.ifEmpty([])
    file ('alignment/sample/picard_metrics/*') from merge_sample_metrics_mqc.collect().ifEmpty([])
    file ('macs/sample/*') from sample_macs_peak_mqc.collect().ifEmpty([])
    file ('macs/sample/*') from sample_macs_qc_mqc.collect().ifEmpty([])
    file ('macs/sample/consensus/*') from sample_macs_consensus_counts_mqc.collect().ifEmpty([])
    file ('macs/sample/consensus/*') from sample_macs_consensus_deseq_mqc.collect().ifEmpty([])
    file ('software_versions/*') from software_versions_yaml.collect()
    file ('workflow_summary/*') from create_workflow_summary(summary)

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc . -f $rtitle $rfilename --config $multiqc_config \\
        -m custom_content -m fastqc -m cutadapt -m samtools -m picard -m featureCounts
    """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                             IGV                                     -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * STEP 8 - Create IGV session file
 */
process igv {
    publishDir "${params.outdir}/igv", mode: 'copy'

    input:
    file ('bwa/replicate/bigwig/*') from replicate_bigwig_igv.collect()
    file ('bwa/replicate/macs/*') from replicate_macs_peaks_igv.collect{it[1]}.ifEmpty([])
    file ('bwa/replicate/macs/consensus/*') from replicate_macs_consensus_igv.collect().ifEmpty([])
    file ('bwa/replicate/macs/consensus/deseq2/*') from replicate_macs_consensus_deseq_comp_bed.collect().ifEmpty([])

    file ('bwa/sample/bigwig/*') from sample_bigwig_igv.collect().ifEmpty([])
    file ('bwa/sample/macs/*') from sample_macs_peaks_igv.collect{it[1]}.ifEmpty([])
    file ('bwa/sample/macs/consensus/*') from sample_macs_consensus_igv.collect().ifEmpty([])
    file ('bwa/sample/macs/consensus/deseq2/*') from sample_macs_consensus_deseq_comp_bed.collect().ifEmpty([])

    output:
    file "*.xml" into igv_session

    script: // scripts are bundled with the pipeline, in nf-core/atacseq/bin/
    """
    igv_get_files.sh ./ mSm > igv_files.txt
    igv_get_files.sh ./ mRp >> igv_files.txt
    igv_files_to_session.py igv_session.xml igv_files.txt ${params.fasta}
    """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       REPORTS/DOCUMENTATION                         -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * STEP 9 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs from output_docs_ch

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/atacseq] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/atacseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/atacseq] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/atacseq] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[nf-core/atacseq] Pipeline Complete"

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                        END OF PIPELINE                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
