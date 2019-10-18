# nf-core/atacseq: Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unpublished Version / DEV]

### `Added`

* [#35](https://github.com/nf-core/atacseq/issues/35) - Add deepTools plotFingerprint
* [#46](https://github.com/nf-core/atacseq/issues/46) - Missing gene_bed path in igenomes config
* Merged in TEMPLATE branch for automated syncing
* Renamed all channels to start with `ch_` prefix
* Capitalised process names
* Add quick start information to main README
* Add parameters:
  * `--seq_center`
  * `--fingerprint_bins`
  * `--broad_cutoff`
  * `--min_reps_consensus`
  * `--save_macs_pileup`
  * `--skip_diff_analysis`
  * `--skip*` for skipping QC steps
* Change parameter `saveGenomeIndex` to `save_reference`
* Update template to tools `1.7`
* Bump Nextflow version to `19.04.0`
* Change all parameters from `camelCase` to `snake_case`
* Change `--design` parameter to `--input` for standardisation

### `Fixed`

* [#41](https://github.com/nf-core/atacseq/issues/41) - Docs: Add example plot images
* [#44](https://github.com/nf-core/atacseq/issues/44) - Output directory missing: macs2/consensus/deseq2
* [#45](https://github.com/nf-core/atacseq/issues/45) - Wrong x-axis scale for the HOMER: Peak annotation Counts tab plot?
* [#46](https://github.com/nf-core/atacseq/issues/46) - Stage blacklist file in channel properly
* [#50](https://github.com/nf-core/atacseq/issues/50) - HOMER number of peaks does not correspond to found MACS2 peaks
* Smaller logo for completion email
* IGV uses relative instead of absolute paths

### `Dependencies`

* Added preseq v2.0.3
* Added deeptools v3.2.1
* Added r-xfun v0.3
* Added gawk v4.2.1

## [1.0.0] - 2019-04-09

### `Added`
Initial release of nf-core/atacseq
Created with version 1.1 of the [nf-core](http://nf-co.re/) template.
