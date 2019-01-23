# nf-core/atacseq: Usage

## Table of contents

* [Introduction](#general-nextflow-info)
* [Running the pipeline](#running-the-pipeline)
* [Updating the pipeline](#updating-the-pipeline)
* [Reproducibility](#reproducibility)
* [Main arguments](#main-arguments)
    * [`-profile`](#-profile-single-dash)
        * [`awsbatch`](#awsbatch)
        * [`conda`](#conda)
        * [`docker`](#docker)
        * [`singularity`](#singularity)
        * [`test`](#test)
    * [`--design`](#--design)
* [Generic arguments](#generic-arguments)
    * [`--singleEnd`](#--singleEnd)
    * [`--narrowPeak`](#--narrowPeak)
    * [`--fragment_size`](#--fragment_size)
* [Reference genomes](#reference-genomes)
    * [`--genome`](#--genome)
    * [`--fasta`](#--fasta)
    * [`--gtf`](#--gtf)
    * [`--bwa_index_dir`](#--bwa_index_dir)
    * [`--bwa_index_base`](#--bwa_index_base)
    * [`--gene_bed`](#--gene_bed)
    * [`--tss_bed`](#--tss_bed)
    * [`--mito_name`](#--mito_name)
    * [`--macs_gsize`](#--macs_gsize)
    * [`--blacklist`](#--blacklist)
    * [`--saveReference`](#--saveReference)
    * [`--igenomesIgnore`](#--igenomesignore)
* [Adapter trimming](#adapter-trimming)
    * [`--skipTrimming`](#--skipTrimming)
    * [`--saveTrimmed`](#--saveTrimmed)
* [Alignments](#alignments)
    * [`--keepMito`](#--keepMito)
    * [`--keepDups`](#--keepDups)
    * [`--keepMultiMap`](#--keepMultiMap)
    * [`--skipMergeBySample`](#--skipMergeBySample)
    * [`--saveAlignedIntermediates`](#--saveAlignedIntermediates)
* [Job resources](#job-resources)
* [Automatic resubmission](#automatic-resubmission)
* [Custom resource requests](#custom-resource-requests)
* [AWS batch specific parameters](#aws-batch-specific-parameters)
    * [`-awsbatch`](#-awsbatch)
    * [`--awsqueue`](#--awsqueue)
    * [`--awsregion`](#--awsregion)
* [Other command line parameters](#other-command-line-parameters)
    * [`--outdir`](#--outdir)
    * [`--email`](#--email)
    * [`-name`](#-name-single-dash)
    * [`-resume`](#-resume-single-dash)
    * [`-c`](#-c-single-dash)
    * [`--custom_config_version`](#--custom_config_version)
    * [`--max_memory`](#--max_memory)
    * [`--max_time`](#--max_time)
    * [`--max_cpus`](#--max_cpus)
    * [`--plaintext_email`](#--plaintext_email)
    * [`--multiqc_config`](#--multiqc_config)

## General Nextflow info
Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline
The typical command for running the pipeline is as follows:
```bash
nextflow run nf-core/atacseq --design design.csv --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/atacseq
```

### Reproducibility
It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/atacseq releases page](https://github.com/nf-core/atacseq/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Main arguments

### `-profile`
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile docker` - the order of arguments is important!

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

* `awsbatch`
    * A generic configuration profile to be used with AWS Batch.
* `conda`
    * A generic configuration profile to be used with [conda](https://conda.io/docs/)
    * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `docker`
    * A generic configuration profile to be used with [Docker](http://docker.com/)
    * Pulls software from dockerhub: [`nfcore/atacseq`](http://hub.docker.com/r/nfcore/atacseq/)
* `singularity`
    * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
    * Pulls software from singularity-hub
* `test`
    * A profile with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters

### `--design`
You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location.

```bash
--design '[path to design file]'
```

It has to be a comma-separated file with 4 columns, and a header row as shown in the example below:

```bash
sample,replicate,fastq_1,fastq_2
control,1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
control,2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
treatment,1,AEG588A3_S5_L003_R1_001.fastq.gz,AEG588A3_S5_L003_R2_001.fastq.gz
treatment,2,AEG588A4_S3_L002_R1_001.fastq.gz,AEG588A4_S3_L002_R2_001.fastq.gz
treatment,2,AEG588A5_S4_L002_R1_001.fastq.gz,AEG588A5_S4_L002_R2_001.fastq.gz
```

| Column      | Description                                                                                   |
|-------------|-----------------------------------------------------------------------------------------------|
| `sample`    | Group identifier for sample.                                                                  |
| `replicate` | Integer representing replicate number. Must start from 1..<number of replicates>              |
| `fastq_1`   | Full path to FastQ file for read 1. File has to be zipped and have the extension ".fastq.gz". |
| `fastq_2`   | Full path to FastQ file for read 2. File has to be zipped and have the extension ".fastq.gz". |

If you have sequenced the same library more than once you just provide this as a separate entry in the design file with the same replicate identifier. The alignments will be performed separately, and subsequently merged before further analysis.

Example design files have been provided with the pipeline for [paired-end](../assets/design_pe.csv) and [single-end](../assets/design_se.csv) data.

>NB: The pipeline will assume that the design file is in the correct format so please double-check this before execution.

## Generic arguments

### `--singleEnd`
By default, the pipeline expects paired-end data. If you have single-end data, specify `--singleEnd` on the command line when you launch the pipeline.

It is not possible to run a mixture of single-end and paired-end files in one run.

### `--narrowPeak`
MACS2 is run by default with the [`--broad`](https://github.com/taoliu/MACS#--broad) flag. Specify this flag to call peaks in narrowPeak mode.

### `--fragment_size`
Number of base pairs to extend single-end reads when creating bigWig files. Default: `0`

## Reference genomes

The pipeline config files come bundled with paths to the illumina iGenomes reference index files. If running with docker or AWS, the configuration is set up to use the [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) resource.

### `--genome` (using iGenomes)
There are 31 different species supported in the iGenomes references. To run the pipeline, you must specify which to use with the `--genome` flag.

You can find the keys to specify the genomes in the [`iGenomes config file`](../conf/igenomes.config). Common genomes that are supported are:

* Human
  * `--genome GRCh37`
* Mouse
  * `--genome GRCm38`
* _Drosophila_
  * `--genome BDGP6`
* _S. cerevisiae_
  * `--genome 'R64-1-1'`

> There are numerous others - check the config file for more.

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the iGenomes resource. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'GRCh37' {
      fasta   = '<path to the genome fasta file>' // Used if no star index given
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

### `--fasta`
Full path to fasta file containing reference genome (*mandatory* if `--genome` is not specified). If you don't have a BWA index available this will be generated for you automatically. Combine with `--saveReference` to save BWA index for future runs.
```bash
--fasta '[path to FASTA reference]'
```

### `--gtf`
The full path to GTF file for annotating peaks (*mandatory* if `--genome` is not specified). Note that the GTF file should resemble the Ensembl format.
```bash
--gtf '[path to GTF file]'
```

### `--bwa_index_dir`
Directory containing an existing BWA index for your reference genome.
```bash
--bwa_index_dir '[directory containing BWA index]'
```

### `--bwa_index_base`
Base file name for an existing BWA index for your reference genome. Default: `genome.fa`.
```bash
--bwa_index_base '[basename of BWA index]'
```

### `--gene_bed`
The full path to BED file for genome-wide gene intervals. This will be created from the GTF file if it isnt specified.
```bash
--gene_bed '[path to gene BED file]'
```

### `--tss_bed`
The full path to BED file for genome-wide transcription start sites. This will be created from the gene BED file if it isnt specified.
```bash
--tss_bed '[path to tss BED file]'
```

### `--macs_gsize`
[Effective genome size](https://github.com/taoliu/MACS#-g--gsize) parameter required by MACS2. These have been provided when `--genome` is set as *GRCh37*, *GRCm38*, *WBcel235*, *BDGP6*, *R64-1-1*, *EF2*, *hg38*, *hg19* and *mm10*. For other genomes, if this parameter isnt specified then the MACS2 peak-calling and differential analysis will be skipped.
```bash
--macs_gsize 2.7e9
```

### `--mito_name`
Name of Mitochondrial chomosome in genome fasta. Reads aligning to this contig are filtered out if a valid identifier is provided otherwise this step is skipped. Where possible these have been provided in the [`igenomes.config`](../conf/igenomes.config).
```bash
--mito_name chrM
```

### `--blacklist`
If provided, alignments that overlap with the regions in this file will be filtered out (see [ENCODE blacklists](https://sites.google.com/site/anshulkundaje/projects/blacklists)). The file should be in BED format. Blacklisted regions for *GRCh37*, *GRCm38*, *hg19*, *hg38*, *mm10* are bundled with the pipeline in the [`blacklists`](../blacklists/) directory, and as such will be automatically used if any of those genomes are specified with the `--genome` parameter.
```bash
--blacklist '[path to blacklisted regions]'
```

### `--saveReference`
Supply this parameter to save any generated reference genome files such as the BWA index to your results folder. These can then be used for future pipeline runs, reducing processing times.

### `--igenomesIgnore`
Do not load `igenomes.config` when running the pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.

## Adapter trimming
The pipeline accepts a number of parameters to change how the trimming is done, according to your data type.
You can specify custom trimming parameters as follows:

* `--clip_r1 <NUMBER>`
  * Instructs Trim Galore to remove bp from the 5' end of read 1 (for single-end reads).
* `--clip_r2 <NUMBER>`
  * Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only).
* `--three_prime_clip_r1 <NUMBER>`
  * Instructs Trim Galore to remove bp from the 3' end of read 1 _AFTER_ adapter/quality trimming has been
* `--three_prime_clip_r2 <NUMBER>`
  * Instructs Trim Galore to re move bp from the 3' end of read 2 _AFTER_ adapter/quality trimming has been performed.

### `--skipTrimming`
Skip the adapter trimming step. Use this if your input FastQ files have already been trimmed outside of the workflow or if you're very confident that there is no adapter contamination in your data.

### `--saveTrimmed`
By default, trimmed FastQ files will not be saved to the results directory. Specify this flag (or set to true in your config file) to copy these files to the results directory when complete.

## Alignments

### `--keepMito`
Reads mapping to mitochondrial contig are not filtered from alignments.

### `--keepDups`
Duplicate reads are not filtered from alignments.

### `--keepMultiMap`
Reads mapping to multiple locations in the genome are not filtered from alignments.

### `--skipMergeBySample`
Do not perform alignment merging and downstream analysis at the sample-level i.e. only do this at the replicate-level.

### `--saveAlignedIntermediates`
By default, intermediate BAM files will not be saved. The final BAM files created after the appropriate filtering step are always saved to limit storage usage. Set to true to also copy out BAM files from BWA and sorting/filtering steps.

## Job resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests
Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [`Slack`](https://nf-core-invite.herokuapp.com/).

## AWS Batch specific parameters
Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use the `-awsbatch` profile and then specify all of the following parameters.
### `--awsqueue`
The JobQueue that you intend to use on AWS Batch.
### `--awsregion`
The AWS region to run your job in. Default is set to `eu-west-1` but can be adjusted to your needs.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

### `--outdir`
The output directory where the results will be saved.

### `--email`
Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `-name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`
Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`
Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default is set to `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--max_memory`
Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`
Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`
Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`
Set to receive plain-text e-mails instead of HTML formatted.

### `--multiqc_config`
Specify a path to a custom MultiQC configuration file.
