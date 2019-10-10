# nf-core/atacseq: Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unpublished Version / DEV]

### `Added`

* [#35](https://github.com/nf-core/atacseq/issues/35) - Add deepTools plotFingerprint
* [#46](https://github.com/nf-core/atacseq/issues/46) - Missing gene_bed path in igenomes config
* [#46](https://github.com/nf-core/atacseq/issues/46) - Stage blacklist file in channel properly
* Merged in TEMPLATE branch for automated syncing
* Renamed all channels to start with `ch_` prefix
* Capitalised process names
* Add quick start information to main README
* Add parameters:
  * `--seq_center`
  * `--fingerprint_bins`
  * `--broad_cutoff`
  * `--min_reps_consensus`
  * `--saveMACSPileup`
  * `--skipDiffAnalysis`
  * `--skip*` for skipping QC steps
* Update template to tools v1.7

### `Fixed`

* [#41](https://github.com/nf-core/atacseq/issues/41) - Docs: Add example plot images
* [#44](https://github.com/nf-core/atacseq/issues/44) - Output directory missing: macs2/consensus/deseq2
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
