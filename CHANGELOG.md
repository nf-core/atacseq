# nf-core/atacseq: Changelog

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unpublished Version / DEV]

### `Added`

* [#63](https://github.com/nf-core/atacseq/issues/63) - Added multicore support for Trim Galore!
* [#75](https://github.com/nf-core/atacseq/issues/75) - Include gene annotation versions in multiqc report
* [#76](https://github.com/nf-core/atacseq/issues/76) - featureCounts coupled to DESeq2
* [79](https://github.com/nf-core/atacseq/issues/79) - Parallelize DESeq2
* [#80](https://github.com/nf-core/atacseq/pull/80) - Added social preview image
* Update template to tools `1.9`
* Parameter `--skip_consensus_peaks` to skip consensus peak generation

### `Fixed`

* [nf-core/chipseq#118](https://github.com/nf-core/chipseq/issues/118) - Running on with SGE
* [nf-core/chipseq#132](https://github.com/nf-core/chipseq/issues/132) - BigWig Error: sort: cannot create temporary file in '': Read-only file system
* [#73](https://github.com/nf-core/atacseq/issues/73) - macs_annotatePeaks.mLb.clN.summary.txt file is not created
* Make executables in `bin/` compatible with Python 3

### `Dependencies`

* Add python `3.7.6`
* Add markdown `3.1.1`
* Add pymdown-extensions `6.0`
* Add pygments `2.5.2`
* Add pigz `2.3.4`
* Add r-tidyr `1.0.2`
* Add bioconductor-biocparallel `1.20.0`
* Update gawk `4.2.1` -> `5.0.1`
* Update r-base `3.4.1` -> `3.6.2`
* Update r-optparse `1.6.0` -> `1.6.4`
* Update r-ggplot2 `3.1.0` -> `3.2.1`
* Update r-pheatmap `1.0.10` -> `1.0.12`
* Update r-lattice `0.20_35` -> `0.20_40`
* Update r-upsetr `1.3.3` -> `1.4.0`
* Update r-scales `1.0.0` -> `1.1.0`  
* Update r-xfun `0.3` -> `0.12`
* Update fastqc `0.11.8` -> `0.11.9`
* Update trim-galore `0.5.0` -> `0.6.5`
* Update picard `2.19.0` -> `2.21.9`
* Update pysam `0.15.2` -> `0.15.4`
* Update bedtools `2.27.1` -> `2.29.2`
* Update ucsc-bedgraphtobigwig `377` -> `357`
* Update deeptools `3.2.1` -> `3.3.2`
* Update macs2 `2.1.2` -> `2.2.6`
* Update homer `4.9.1` -> `4.10`
* Update ataqv `1.0.0` -> `1.1.1`
* Update subread `1.6.4` -> `2.0.0`
* Update multiqc `1.7` -> `1.8`
* Update bioconductor-deseq2 `1.20.0` -> `1.26.0`
* Update bioconductor-vsn `3.46.0` -> `3.54.0`
* Remove r-reshape2 `1.4.3`

### `Deprecated`

## [1.1.0] - 2019-11-05

### `Added`

* [#35](https://github.com/nf-core/atacseq/issues/35) - Add deepTools plotFingerprint
* [#46](https://github.com/nf-core/atacseq/issues/46) - Missing gene_bed path in igenomes config
* Merged in TEMPLATE branch for automated syncing
* Update template to tools `1.7`
* Add `CITATIONS.md` file
* Capitalised process names
* Add parameters:
  * `--seq_center`
  * `--trim_nextseq`
  * `--fingerprint_bins`
  * `--broad_cutoff`
  * `--min_reps_consensus`
  * `--save_macs_pileup`
  * `--skip_diff_analysis`
  * `--skip_*` for skipping QC steps

### `Fixed`

* **Change all parameters from `camelCase` to `snake_case` (see [Deprecated](#Deprecated))**
* [#41](https://github.com/nf-core/atacseq/issues/41) - Docs: Add example plot images
* [#44](https://github.com/nf-core/atacseq/issues/44) - Output directory missing: macs2/consensus/deseq2
* [#45](https://github.com/nf-core/atacseq/issues/45) - Wrong x-axis scale for the HOMER: Peak annotation Counts tab plot?
* [#46](https://github.com/nf-core/atacseq/issues/46) - Stage blacklist file in channel properly
* [#50](https://github.com/nf-core/atacseq/issues/50) - HOMER number of peaks does not correspond to found MACS2 peaks
* Fixed bug in UpSetR peak intersection plot
* IGV now uses relative instead of absolute paths
* Smaller logo for completion email
* Renamed all channels to start with `ch_` prefix
* Increase default resource requirements in `base.config`
* Increase process-specific requirements based on user-reported failures

### `Dependencies`

* Update Nextflow `0.32.0` -> `19.10.0`
* Add preseq `2.0.3`
* Add deeptools `3.2.1`
* Add r-xfun `0.3`
* Add gawk `4.2.1`

### `Deprecated`

| Deprecated                   | Replacement               |
|------------------------------|---------------------------|
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
