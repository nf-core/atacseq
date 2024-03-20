# nf-core/atacseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.2.0dev - [date]

### Enhancements & fixes

- [[#327](https://github.com/nf-core/atacseq/issues/327)]Consistently support `.csi` indices as alternative to `.bai` to allow SAMTOOLS_INDEX to be used with the `-c` flag.
- Updated pipeline template to [nf-core/tools 2.10](https://github.com/nf-core/tools/releases/tag/2.10)

## [[2.1.2](https://github.com/nf-core/atacseq/releases/tag/2.1.2)] - 2022-08-07

### Enhancements & fixes

- [[#322](https://github.com/nf-core/atacseq/issues/322)]Remove fasta from required schema parameters so that when launching from tools it is not required.
- Updates `homer_annotatepeaks` module.

## [[2.1.1](https://github.com/nf-core/atacseq/releases/tag/2.1.1)] - 2022-07-21

- Minor patch release to fix AWS full test.

## [[2.1.0](https://github.com/nf-core/atacseq/releases/tag/2.1.0)] - 2022-07-21

### Credits

Special thanks to the following for their contributions to the release:

- [Adam Talbot](https://github.com/adamrtalbot)
- [Björn Langer](https://github.com/bjlang)
- [Harshil Patel](https://github.com/drpatelh)
- [Maxime Garcia](https://github.com/maxulysse)
- [Rob Syme](https://github.com/robsyme)

Thank you to everyone else that has contributed by reporting bugs, enhancements or in any other way, shape or form.

### Enhancements & fixes

- [[#226](https://github.com/nf-core/atacseq/issues/262)] - Add `ataqv_mito_reference` parameter.
- Optional support of control data analog to nf-core/chipseq.
- [[#277](https://github.com/nf-core/atacseq/issues/277)] - Fix error when using a gunziped fasta file.
- [[#286](https://github.com/nf-core/atacseq/issues/286)] - Fix error when no `--mito_name parameter is provided.
- [[#268](https://github.com/nf-core/atacseq/issues/268)] - Fix error when a bed file is provided using the `--blacklist` option.
- [[#278](https://github.com/nf-core/atacseq/issues/278)] - Make genome fasta file available when `IGV` process is run.
- [[#276](https://github.com/nf-core/atacseq/issues/276)] - Bump version 1.3.1 of ataqv to fix enrichment plots rendering.
- [[#290](https://github.com/nf-core/atacseq/issues/290)] - Fix case-sensitivity issue while sorting bedGraph.
- [[#295](https://github.com/nf-core/atacseq/issues/295)] - Enable downstream steps for bam files produced from paired-end reads by `chromap` after its upgrade.
- Updated pipeline template to [nf-core/tools 2.9](https://github.com/nf-core/tools/releases/tag/2.9)
- Make fasta index available for IGV session.
- [[nf-core/chipseq#347](https://github.com/nf-core/chipseq/issues/347)] - Add read group tag to bam files processed by bowtie2

### Parameters

| Old parameter | New parameter            |
| ------------- | ------------------------ |
|               | `--with_control`         |
|               | `--ataqv_mito_reference` |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
> **NB:** Parameter has been **added** if just the new parameter information is present.
> **NB:** Parameter has been **removed** if parameter information isn't present.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency              | Old version | New version |
| ----------------------- | ----------- | ----------- |
| `ataqv`                 | 1.3.0       | 1.3.1       |
| `chromap`               | 0.2.1       | 0.2.4       |
| `multiqc`               | 1.13        | 1.14        |
| `picard`                | 2.27.4      | 3.0.0       |
| `samtools`              | 1.15.1      | 1.17        |
| `ucsc-bedgraphtobigwig` | 377         | 445         |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

## [[2.0](https://github.com/nf-core/atacseq/releases/tag/2.0)] - 2022-11-30

### :warning: Major enhancements

- Samplesheet format has changed from `group,replicate,fastq_1,fastq_2` to `sample,fastq_1,fastq_2,replicate`.
  - This is primarily because we have removed the differential accessibility from the pipeline in favour of using more formal statistical methods downstream with the pipeline results. We kept replicate information as mandatory to merge samples at the replicate level to improve coverage for footprinting analysis.
- Pipeline has been re-implemented in [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html)
- Bump minimum Nextflow version from `19.10.0` -> `21.10.3`
- Updated pipeline template to [nf-core/tools 2.6](https://github.com/nf-core/tools/releases/tag/2.6)
- [[#135](https://github.com/nf-core/atacseq/issues/135)] - ERROR: Please check design file header
- [[#158](https://github.com/nf-core/atacseq/issues/158)] - Problem using hyphens in group names
- [[#182](https://github.com/nf-core/atacseq/issues/182)] - Update `macs_gsize` in `igenomes.config`, create a new `--read_length` parameter and implement the logic to calculate `--macs_gsize` when the parameter is missing
- [[#193](https://github.com/nf-core/atacseq/issues/193)] - TSS heatmap inconsistency between ChIPseeker and nf-core ATAC-seq pipeline
- [[#201](https://github.com/nf-core/atacseq/issues/201)] - Update blacklist bed files.
- [[#220](https://github.com/nf-core/atacseq/issues/220)] - Move FastQ files for minimal tests to S3
- [[nf-core/chipseq#160](https://github.com/nf-core/chipseq/issues/160)] - Add `bowtie2` and `star` as available aligners, via the `--aligner` parameter
- [[nf-core/chipseq#233](https://github.com/nf-core/chipseq/issues/233)] - Add `chromap` to the available aligners
- Add `--save_unaligned` parameter (only available for `bowtie2` and `star`)
- Preseq will be skipped by default (see [nf-core/rnaseq#897](https://github.com/nf-core/rnaseq/issues/897))
- `--deseq2_vst` will be turned on by default

### Parameters

| Old parameter          | New parameter      |
| ---------------------- | ------------------ |
| `--clusterOptions`     |                    |
| `--conda`              | `--enable_conda`   |
|                        | `--skip_qc`        |
|                        | `--aligner`        |
|                        | `--save_unaligned` |
|                        | `--bowtie2_index`  |
|                        | `--chromap_index`  |
|                        | `--star_index`     |
| `--skip_diff_analysis` | `--skip_deseq2_qc` |
| `--single_end`         |                    |
|                        | `--read_length`    |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
> **NB:** Parameter has been **added** if just the new parameter information is present.
> **NB:** Parameter has been **removed** if parameter information isn't present.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency              | Old version | New version |
| ----------------------- | ----------- | ----------- |
| `ataqv`                 | 1.1.1       | 1.3.0       |
| `bamtools`              | 2.5.1       | 2.5.2       |
| `bedtools`              | 2.29.2      | 2.30.0      |
| `bioconductor-deseq2`   | 1.26.0      | 1.28.0      |
| `deeptools`             | 3.4.3       | 3.5.1       |
| `multiqc`               | 1.9         | 1.13        |
| `ucsc-bedgraphtobigwig` | 357         | 377         |
| `samtools`              | 1.10        | 1.16.1      |
| `picard`                | 2.23.1      | 2.27.4      |
| `preseq`                | 2.0.3       | 3.1.2       |
| `pysam`                 | 0.15.3      | 0.19.0      |
| `r-base`                | 3.6.1       | 4.0.3       |
| `r-ggplot2`             | 3.3.2       | 3.3.3       |
| `trim-galore`           | 0.6.5       | 0.6.7       |
| `khmer`                 | -           | 3.0.0a3     |
| `r-optparse`            | -           | 1.7.1       |
| `r-tidyr`               | -           | -           |
| `r-lattice`             | -           | -           |
| `r-xfun`                | -           | -           |
| `bioconductor-vsn`      | -           | -           |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

## [1.2.2] - 2022-05-12

- Minor patch release to fix Conda environment and enable DSL1 by default in the `nextflow.config`.

### `Dependencies`

- Update r-base `3.6.2` -> `3.6.3`
- Update r-xfun `0.15` -> `0.22`

## [1.2.1] - 2020-07-29

- [#118](https://github.com/nf-core/atacseq/issues/118) - Minor patch release to update pipeline schema

## [1.2.0] - 2020-07-02

### `Added`

- [#63](https://github.com/nf-core/atacseq/issues/63) - Added multicore support for Trim Galore!
- [#75](https://github.com/nf-core/atacseq/issues/75) - Include gene annotation versions in multiqc report
- [#76](https://github.com/nf-core/atacseq/issues/76) - featureCounts coupled to DESeq2
- [#79](https://github.com/nf-core/atacseq/issues/79) - Parallelize DESeq2
- [#80](https://github.com/nf-core/atacseq/pull/80) - Added social preview image
- [#97](https://github.com/nf-core/atacseq/issues/97) - PBC1, PBC2 from pipeline?
- [#107](https://github.com/nf-core/atacseq/issues/107) - Add options to change MACS2 parameters
- [nf-core/chipseq#153](https://github.com/nf-core/chipseq/issues/153) - Add plotHeatmap
- [nf-core/chipseq#159](https://github.com/nf-core/chipseq/issues/159) - expose bwa mem -T parameter
- Regenerated screenshots and added collapsible sections for output files in `docs/output.md`
- Update template to tools `1.9`
- Replace `set` with `tuple` and `file()` with `path()` in all processes
- Capitalise process names
- Parameters:
  - `--bwa_min_score` to set minimum alignment score for BWA MEM
  - `--macs_fdr` to provide FDR threshold for MACS2 peak calling
  - `--macs_pvalue` to provide p-value threshold for MACS2 peak calling
  - `--skip_peak_qc` to skip MACS2 peak QC plot generation
  - `--skip_peak_annotation` to skip annotation of MACS2 and consensus peaks with HOMER
  - `--skip_consensus_peaks` to skip consensus peak generation
  - `--deseq2_vst` to use variance stabilizing transformation (VST) instead of regularized log transformation (rlog) with DESeq2
  - `--publish_dir_mode` to customise method of publishing results to output directory [nf-core/tools#585](https://github.com/nf-core/tools/issues/585)

### `Fixed`

- [#71](https://github.com/nf-core/atacseq/issues/71) - consensus_peaks.mLb.clN.boolean.intersect.plot.pdf not generated
- [#73](https://github.com/nf-core/atacseq/issues/73) - macs_annotatePeaks.mLb.clN.summary.txt file is not created
- [#86](https://github.com/nf-core/atacseq/issues/86) - bug in the plot_homer_annotatepeaks.r script
- [#102](https://github.com/nf-core/atacseq/issues/102) - Incorrect Group ID assigned by featurecounts_deseq2.r
- [#110](https://github.com/nf-core/atacseq/pull/110) - updated AWS test GitHub actions
- [#109](https://github.com/nf-core/atacseq/issues/109) - Specify custom gtf but gene bed is not generated from that gtf?
- [nf-core/chipseq#118](https://github.com/nf-core/chipseq/issues/118) - Running on with SGE
- [nf-core/chipseq#132](https://github.com/nf-core/chipseq/issues/132) - BigWig Error: sort: cannot create temporary file in '': Read-only file system
- [nf-core/chipseq#154](https://github.com/nf-core/chipseq/issues/154) - computeMatrix.val.mat.gz files not zipped
- Make executables in `bin/` compatible with Python 3

### `Dependencies`

- Add bioconductor-biocparallel `1.20.0`
- Add markdown `3.2.2`
- Add pigz `2.3.4`
- Add pygments `2.6.1`
- Add pymdown-extensions `7.1`
- Add python `3.7.6`
- Add r-reshape2 `1.4.4`
- Add r-tidyr `1.1.0`
- Update ataqv `1.0.0` -> `1.1.1`
- Update bedtools `2.27.1` -> `2.29.2`
- Update bioconductor-deseq2 `1.20.0` -> `1.26.0`
- Update bioconductor-vsn `3.46.0` -> `3.54.0`
- Update deeptools `3.2.1` -> `3.4.3`
- Update fastqc `0.11.8` -> `0.11.9`
- Update gawk `4.2.1` -> `5.1.0`
- Update homer `4.9.1` -> `4.11`
- Update macs2 `2.1.2` -> `2.2.7.1`
- Update multiqc `1.7` -> `1.8`
- Update picard `2.19.0` -> `2.23.1`
- Update pysam `0.15.2` -> `0.15.3`
- Update r-base `3.4.1` -> `3.6.2`
- Update r-ggplot2 `3.1.0` -> `3.3.2`
- Update r-lattice `0.20_35` -> `0.20_41`
- Update r-optparse `1.6.0` -> `1.6.6`
- Update r-pheatmap `1.0.10` -> `1.0.12`
- Update r-scales `1.0.0` -> `1.1.1`
- Update r-upsetr `1.3.3` -> `1.4.0`
- Update r-xfun `0.3` -> `0.15`
- Update samtools `1.9` -> `1.10`
- Update subread `1.6.4` -> `2.0.1`
- Update trim-galore `0.5.0` -> `0.6.5`
- Update ucsc-bedgraphtobigwig `377` -> `357`

## [1.1.0] - 2019-11-05

### `Added`

- [#35](https://github.com/nf-core/atacseq/issues/35) - Add deepTools plotFingerprint
- [#46](https://github.com/nf-core/atacseq/issues/46) - Missing gene_bed path in igenomes config
- Merged in TEMPLATE branch for automated syncing
- Update template to tools `1.7`
- Add `CITATIONS.md` file
- Capitalised process names
- Add parameters:
  - `--seq_center`
  - `--trim_nextseq`
  - `--fingerprint_bins`
  - `--broad_cutoff`
  - `--min_reps_consensus`
  - `--save_macs_pileup`
  - `--skip_diff_analysis`
  - `--skip_*` for skipping QC steps

### `Fixed`

- **Change all parameters from `camelCase` to `snake_case` (see [Deprecated](#deprecated))**
- [#41](https://github.com/nf-core/atacseq/issues/41) - Docs: Add example plot images
- [#44](https://github.com/nf-core/atacseq/issues/44) - Output directory missing: macs2/consensus/deseq2
- [#45](https://github.com/nf-core/atacseq/issues/45) - Wrong x-axis scale for the HOMER: Peak annotation Counts tab plot?
- [#46](https://github.com/nf-core/atacseq/issues/46) - Stage blacklist file in channel properly
- [#50](https://github.com/nf-core/atacseq/issues/50) - HOMER number of peaks does not correspond to found MACS2 peaks
- Fixed bug in UpSetR peak intersection plot
- IGV now uses relative instead of absolute paths
- Smaller logo for completion email
- Renamed all channels to start with `ch_` prefix
- Increase default resource requirements in `base.config`
- Increase process-specific requirements based on user-reported failures

### `Dependencies`

- Update Nextflow `0.32.0` -> `19.10.0`
- Add preseq `2.0.3`
- Add deeptools `3.2.1`
- Add r-xfun `0.3`
- Add gawk `4.2.1`

### `Deprecated`

| Deprecated                   | Replacement               |
| ---------------------------- | ------------------------- |
| `--design`                   | `--input`                 |
| `--singleEnd`                | `--single_end`            |
| `--saveGenomeIndex`          | `--save_reference`        |
| `--skipTrimming`             | `--skip_trimming`         |
| `--saveTrimmed`              | `--save_trimmed`          |
| `--keepMito`                 | `--keep_mito`             |
| `--keepDups`                 | `--keep_dups`             |
| `--keepMultiMap`             | `--keep_multi_map`        |
| `--skipMergeReplicates`      | `--skip_merge_replicates` |
| `--saveAlignedIntermediates` | `--save_align_intermeds`  |
| `--narrowPeak`               | `--narrow_peak`           |

## [1.0.0] - 2019-04-09

### `Added`

Initial release of nf-core/atacseq
Created with version 1.1 of the [nf-core](http://nf-co.re/) template.
