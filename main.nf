#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/atacseq
========================================================================================
    Github : https://github.com/nf-core/atacseq
    Website: https://nf-co.re/atacseq
    Slack  : https://nfcore.slack.com/channels/atacseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta      = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.bwa_index  = WorkflowMain.getGenomeAttribute(params, 'bwa')
params.gtf        = WorkflowMain.getGenomeAttribute(params, 'gtf')
params.gff        = WorkflowMain.getGenomeAttribute(params, 'gff')
params.gene_bed   = WorkflowMain.getGenomeAttribute(params, 'gene_bed')
params.macs_gsize = WorkflowMain.getGenomeAttribute(params, 'macs_gsize')
params.blacklist  = WorkflowMain.getGenomeAttribute(params, 'blacklist')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { ATACSEQ } from './workflows/atacseq'

//
// WORKFLOW: Run main nf-core/atacseq analysis pipeline
//
workflow NFCORE_ATACSEQ {
    ATACSEQ ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_ATACSEQ ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
