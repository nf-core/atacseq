# nf-core/atacseq: Usage

<<<<<<< HEAD
## Table of contents

* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
* [Running the pipeline](#running-the-pipeline)
    * [Updating the pipeline](#updating-the-pipeline)
    * [Reproducibility](#reproducibility)
* [Main arguments](#main-arguments)
    * [`-profile`](#-profile)
    * [`--input`](#--input)
* [Generic arguments](#generic-arguments)
    * [`--single_end`](#--single_end)
    * [`--seq_center`](#--seq_center)
    * [`--fragment_size`](#--fragment_size)
    * [`--fingerprint_bins`](#--fingerprint_bins)
* [Reference genomes](#reference-genomes)
    * [`--genome` (using iGenomes)](#--genome-using-igenomes)
    * [`--fasta`](#--fasta)
    * [`--gtf`](#--gtf)
    * [`--bwa_index`](#--bwa_index)
    * [`--gene_bed`](#--gene_bed)
    * [`--tss_bed`](#--tss_bed)
    * [`--macs_gsize`](#--macs_gsize)
    * [`--blacklist`](#--blacklist)
    * [`--mito_name`](#--mito_name)
    * [`--save_reference`](#--save_reference)
    * [`--igenomes_ignore`](#--igenomes_ignore)
* [Adapter trimming](#adapter-trimming)
    * [`--skip_trimming`](#--skip_trimming)
    * [`--save_trimmed`](#--save_trimmed)
* [Alignments](#alignments)
    * [`--bwa_min_score`](#--bwa_min_score)
    * [`--keep_mito`](#--keep_mito)
    * [`--keep_dups`](#--keep_dups)
    * [`--keep_multi_map`](#--keep_multi_map)
    * [`--skip_merge_replicates`](#--skip_merge_replicates)
    * [`--save_align_intermeds`](#--save_align_intermeds)
* [Peaks](#peaks)
    * [`--narrow_peak`](#--narrow_peak)
    * [`--broad_cutoff`](#--broad_cutoff)
    * [`--macs_fdr`](#--macs_fdr)
    * [`--macs_pvalue`](#--macs_pvalue)
    * [`--min_reps_consensus`](#--min_reps_consensus)
    * [`--save_macs_pileup`](#--save_macs_pileup)
    * [`--skip_peak_qc`](#--skip_peak_qc)
    * [`--skip_peak_annotation`](#--skip_peak_annotation)
    * [`--skip_consensus_peaks`](#--skip_consensus_peaks)
* [Differential analysis](#differential_analysis)
    * [`--deseq2_vst`](#--deseq2_vst)
    * [`--skip_diff_analysis`](#--skip_diff_analysis)
* [Skipping QC steps](#skipping-qc-steps)
* [Job resources](#job-resources)
    * [Automatic resubmission](#automatic-resubmission)
    * [Custom resource requests](#custom-resource-requests)
* [AWS Batch specific parameters](#aws-batch-specific-parameters)
    * [`--awsqueue`](#--awsqueue)
    * [`--awsregion`](#--awsregion)
    * [`--awscli`](#--awscli)
* [Other command line parameters](#other-command-line-parameters)
    * [`--outdir`](#--outdir)
    * [`--publish_dir_mode`](#--publish_dir_mode)
    * [`--email`](#--email)
    * [`--email_on_fail`](#--email_on_fail)
    * [`--max_multiqc_email_size`](#--max_multiqc_email_size)
    * [`-name`](#-name)
    * [`-resume`](#-resume)
    * [`-c`](#-c)
    * [`--custom_config_version`](#--custom_config_version)
    * [`--custom_config_base`](#--custom_config_base)
    * [`--max_memory`](#--max_memory)
    * [`--max_time`](#--max_time)
    * [`--max_cpus`](#--max_cpus)
    * [`--plaintext_email`](#--plaintext_email)
    * [`--monochrome_logs`](#--monochrome_logs)
    * [`--multiqc_config`](#--multiqc_config)

## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
=======
## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->
>>>>>>> TEMPLATE

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
<<<<<<< HEAD
nextflow run nf-core/atacseq --input design.csv --genome GRCh37 -profile docker
=======
nextflow run nf-core/atacseq --input '*_R{1,2}.fastq.gz' -profile docker
>>>>>>> TEMPLATE
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

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
<<<<<<< HEAD
    * A generic configuration profile to be used with [Docker](http://docker.com/)
    * Pulls software from dockerhub: [`nfcore/atacseq`](http://hub.docker.com/r/nfcore/atacseq/)
* `singularity`
    * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
    * Pulls software from DockerHub: [`nfcore/atacseq`](http://hub.docker.com/r/nfcore/atacseq/)
=======
  * A generic configuration profile to be used with [Docker](https://docker.com/)
  * Pulls software from Docker Hub: [`nfcore/atacseq`](https://hub.docker.com/r/nfcore/atacseq/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
  * Pulls software from Docker Hub: [`nfcore/atacseq`](https://hub.docker.com/r/nfcore/atacseq/)
>>>>>>> TEMPLATE
* `conda`
    * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker or Singularity.
    * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
    * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
    * A profile with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters

### `--input`

<<<<<<< HEAD
You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 4 columns, and a header row as shown in the examples below.

```bash
--input '[path to design file]'
```

#### Multiple replicates

The `group` identifier is the same when you have multiple replicates from the same experimental group, just increment the `replicate` identifier appropriately. The first replicate value for any given experimental group must be 1. Below is an example for a single experimental group in triplicate:

```bash
group,replicate,fastq_1,fastq_2
control,1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
control,2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
control,3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
```

#### Multiple runs of the same library

The `group` and `replicate` identifiers are the same when you have re-sequenced the same sample more than once (e.g. to increase sequencing depth). The pipeline will perform the alignments in parallel, and subsequently merge them before further analysis. Below is an example for two samples sequenced across multiple lanes:

```bash
group,replicate,fastq_1,fastq_2
control,1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
control,1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
treatment,1,AEG588A4_S4_L003_R1_001.fastq.gz,AEG588A4_S4_L003_R2_001.fastq.gz
treatment,1,AEG588A4_S4_L004_R1_001.fastq.gz,AEG588A4_S4_L004_R2_001.fastq.gz
```

#### Full design

A final design file may look something like the one below. This is for two experimental groups in triplicate, where the last replicate of the `treatment` group has been sequenced twice.

```bash
group,replicate,fastq_1,fastq_2
control,1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
control,2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
control,3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
treatment,1,AEG588A4_S4_L003_R1_001.fastq.gz,AEG588A4_S4_L003_R2_001.fastq.gz
treatment,2,AEG588A5_S5_L003_R1_001.fastq.gz,AEG588A5_S5_L003_R2_001.fastq.gz
treatment,3,AEG588A6_S6_L003_R1_001.fastq.gz,AEG588A6_S6_L003_R2_001.fastq.gz
treatment,3,AEG588A6_S6_L004_R1_001.fastq.gz,AEG588A6_S6_L004_R2_001.fastq.gz
```

| Column      | Description                                                                                                 |
|-------------|-------------------------------------------------------------------------------------------------------------|
| `group`     | Group identifier for sample. This will be identical for replicate samples from the same experimental group. |
| `replicate` | Integer representing replicate number. Must start from `1..<number of replicates>`.                         |
| `fastq_1`   | Full path to FastQ file for read 1. File has to be zipped and have the extension ".fastq.gz" or ".fq.gz".   |
| `fastq_2`   | Full path to FastQ file for read 2. File has to be zipped and have the extension ".fastq.gz" or ".fq.gz".   |

Example design files have been provided with the pipeline for [paired-end](../assets/design_pe.csv) and [single-end](../assets/design_se.csv) data.

## Generic arguments

### `--single_end`

By default, the pipeline expects paired-end data. If you have single-end data, specify `--single_end` on the command line when you launch the pipeline.

It is not possible to run a mixture of single-end and paired-end files in one run.

### `--seq_center`

Sequencing center information that will be added to read groups in BAM files.

### `--fragment_size`

Number of base pairs to extend single-end reads when creating bigWig files (Default: `0`).

### `--fingerprint_bins`

Number of genomic bins to use when generating the deepTools fingerprint plot. Larger numbers will give a smoother profile, but take longer to run (Default: `500000`).

## Reference genomes

The pipeline config files come bundled with paths to the illumina iGenomes reference index files. If running with docker or AWS, the configuration is set up to use the [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) resource.

### `--genome` (using iGenomes)
=======
### `-resume`
>>>>>>> TEMPLATE

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

<<<<<<< HEAD
* Human
    * `--genome GRCh37`
* Mouse
    * `--genome GRCm38`
* _Drosophila_
    * `--genome BDGP6`
* _S. cerevisiae_
    * `--genome 'R64-1-1'`
=======
### `-c`
>>>>>>> TEMPLATE

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

#### Custom resource requests

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

<<<<<<< HEAD
=======
Whilst these default requirements will hopefully work for most people with most data, you may find that you want to customise the compute resources that the pipeline requests. You can do this by creating a custom config file. For example, to give the workflow process `star` 32GB of memory, you could use the following config:

>>>>>>> TEMPLATE
```nextflow
process {
  withName: star {
    memory = 32.GB
  }
}
```

<<<<<<< HEAD
### `--fasta`

Full path to fasta file containing reference genome (*mandatory* if `--genome` is not specified). If you don't have a BWA index available this will be generated for you automatically. Combine with `--save_reference` to save BWA index for future runs.

```bash
--fasta '[path to FASTA reference]'
```

### `--gtf`

The full path to GTF file for annotating peaks (*mandatory* if `--genome` is not specified). Note that the GTF file should resemble the Ensembl format.

```bash
--gtf '[path to GTF file]'
```

### `--bwa_index`

Full path to an existing BWA index for your reference genome including the base name for the index.

```bash
--bwa_index '[directory containing BWA index]/genome.fa'
```

### `--gene_bed`

The full path to BED file for genome-wide gene intervals. This will be created from the GTF file if not specified.

```bash
--gene_bed '[path to gene BED file]'
```

### `--tss_bed`

The full path to BED file for genome-wide transcription start sites. This will be created from the gene BED file if not specified.

```bash
--tss_bed '[path to tss BED file]'
```

### `--macs_gsize`

[Effective genome size](https://github.com/taoliu/MACS#-g--gsize) parameter required by MACS2. These have been provided when `--genome` is set as *GRCh37*, *GRCh38*, *GRCm38*, *WBcel235*, *BDGP6*, *R64-1-1*, *EF2*, *hg38*, *hg19* and *mm10*. For other genomes, if this parameter is not specified then the MACS2 peak-calling and differential analysis will be skipped.

```bash
--macs_gsize 2.7e9
```

### `--blacklist`

If provided, alignments that overlap with the regions in this file will be filtered out (see [ENCODE blacklists](https://sites.google.com/site/anshulkundaje/projects/blacklists)). The file should be in BED format. Blacklisted regions for *GRCh37*, *GRCh38*, *GRCm38*, *hg19*, *hg38*, *mm10* are bundled with the pipeline in the [`blacklists`](../assets/blacklists/) directory, and as such will be automatically used if any of those genomes are specified with the `--genome` parameter.

```bash
--blacklist '[path to blacklisted regions]'
```

### `--mito_name`

Name of mitochondrial chomosome in reference assembly. Reads aligning to this contig are filtered out if a valid identifier is provided otherwise this step is skipped. Where possible these have been provided in the [`igenomes.config`](../conf/igenomes.config).

```bash
--mito_name chrM
```

### `--save_reference`

If the BWA index is generated by the pipeline use this parameter to save it to your results folder. These can then be used for future pipeline runs, reducing processing times (Default: false).

### `--igenomes_ignore`

Do not load `igenomes.config` when running the pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config` (Default: false).

## Adapter trimming

The pipeline accepts a number of parameters to change how the trimming is done, according to your data type.
You can specify custom trimming parameters as follows:

* `--clip_r1 [int]`
    * Instructs Trim Galore to remove [int] bp from the 5' end of read 1 (for single-end reads).
* `--clip_r2 [int]`
    * Instructs Trim Galore to remove [int] bp from the 5' end of read 2 (paired-end reads only).
* `--three_prime_clip_r1 [int]`
    * Instructs Trim Galore to remove [int] bp from the 3' end of read 1 _AFTER_ adapter/quality trimming has been
* `--three_prime_clip_r2 [int]`
    * Instructs Trim Galore to remove [int] bp from the 3' end of read 2 _AFTER_ adapter/quality trimming has been performed.
* `--trim_nextseq [int]`
    * This enables the option Cutadapt `--nextseq-trim=3'CUTOFF` option via Trim Galore, which will set a quality cutoff (that is normally given with -q instead), but qualities of G bases are ignored. This trimming is in common for the NextSeq- and NovaSeq-platforms, where basecalls without any signal are called as high-quality G bases.

### `--skip_trimming`

Skip the adapter trimming step. Use this if your input FastQ files have already been trimmed outside of the workflow or if you're very confident that there is no adapter contamination in your data (Default: false).

### `--save_trimmed`

By default, trimmed FastQ files will not be saved to the results directory. Specify this flag (or set to true in your config file) to copy these files to the results directory when complete (Default: false).

## Alignments

### `--bwa_min_score`

Don’t output BWA MEM alignments with score lower than this parameter (Default: false).

### `--keep_mito`

Reads mapping to mitochondrial contig are not filtered from alignments (Default: false).

### `--keep_dups`

Duplicate reads are not filtered from alignments (Default: false).

### `--keep_multi_map`

Reads mapping to multiple locations in the genome are not filtered from alignments (Default: false).

### `--skip_merge_replicates`

An additional series of steps are performed by the pipeline by merging the replicates from the same experimental group. This is primarily to increase the sequencing depth in order to perform downstream analyses such as footprinting. Specifying this parameter means that these steps will not be performed (Default: false).

### `--save_align_intermeds`

By default, intermediate BAM files will not be saved. The final BAM files created after the appropriate filtering step are always saved to limit storage usage. Set to true to also save other intermediate BAM files (Default: false).

## Peaks

### `--narrow_peak`

MACS2 is run by default with the [`--broad`](https://github.com/taoliu/MACS#--broad) flag. Specify this flag to call peaks in narrowPeak mode (Default: false).

### `--broad_cutoff`

Specifies broad cut-off value for MACS2. Only used when `--narrow_peak` isnt specified (Default: `0.1`).

### `--macs_fdr`

Minimum FDR (q-value) cutoff for peak detection, `--macs_fdr` and `--macs_pvalue` are mutually exclusive (Default: false).

### `--macs_pvalue`

p-value cutoff for peak detection, `--macs_fdr` and `--macs_pvalue` are mutually exclusive (Default: false). If `--macs_pvalue` cutoff is set, q-value will not be calculated and reported as -1 in the final .xls file.

### `--min_reps_consensus`

Number of biological replicates required from a given condition for a peak to contribute to a consensus peak . If you are confident you have good reproducibility amongst your replicates then you can increase the value of this parameter to create a "reproducible" set of consensus of peaks. For example, a value of 2 will mean peaks that have been called in at least 2 replicates will contribute to the consensus set of peaks, and as such peaks that are unique to a given replicate will be discarded (Default: 1).

```bash
-- min_reps_consensus 1
```

### `--save_macs_pileup`

Instruct MACS2 to create bedGraph files using the `-B --SPMR` parameters (Default: false).

### `--skip_peak_qc`

Skip MACS2 peak QC plot generation (Default: false).

### `--skip_peak_annotation`

Skip annotation of MACS2 and consensus peaks with HOMER (Default: false).

### `--skip_consensus_peaks`

Skip consensus peak generation, annotation and counting (Default: false).

## Differential analysis

### `--deseq2_vst`

Use `vst` transformation instead of `rlog` with DESeq2. See [DESeq2 docs](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-transformations-and-visualization) (Default: false).

### `--skip_diff_analysis`

Skip differential accessibility analysis with DESeq2 (Default: false).

## Skipping QC steps

The pipeline contains a large number of quality control steps. Sometimes, it may not be desirable to run all of them if time and compute resources are limited.
The following options make this easy:

| Step                      | Description                        |
|---------------------------|------------------------------------|
| `--skip_fastqc`           | Skip FastQC                        |
| `--skip_picard_metrics`   | Skip Picard CollectMultipleMetrics |
| `--skip_preseq`           | Skip Preseq                        |
| `--skip_plot_profile`     | Skip deepTools plotProfile         |
| `--skip_plot_fingerprint` | Skip deepTools plotFingerprint     |
| `--skip_ataqv`            | Skip Ataqv                         |
| `--skip_igv`              | Skip IGV                           |
| `--skip_multiqc`          | Skip MultiQC                       |

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.
=======
See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information.
>>>>>>> TEMPLATE

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

#### Nextflow memory requirements

<<<<<<< HEAD
## Other command line parameters

### `--outdir`

The output directory where the results will be saved.

### `--publish_dir_mode`

Value passed to Nextflow [`publishDir`](https://www.nextflow.io/docs/latest/process.html#publishdir) directive for publishing results in the output directory. Available: 'symlink', 'rellink', 'link', 'copy', 'copyNoFollow' and 'move' (Default: 'copy').

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `--email_on_fail`

This works exactly as with `--email`, except emails are only sent if the workflow is not successful.

### `--max_multiqc_email_size`

Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB).

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

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default: `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:
=======
In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):
>>>>>>> TEMPLATE

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
