# Draft PR: fix module-update regressions after #448

Target: `nf-core/atacseq:dev`  
Related issue: [#473](https://github.com/nf-core/atacseq/issues/473)  
Status: versioned proposal on the fork; bounded module review complete, archive and final PR review pending.

## Problem and changes

The module update in #448 changed scientific output and optional-QC failure behavior. This branch contains four corrections:

- Update the complete Trim Galore module to **2.3.0**, which fixes the upstream wrong-end R2 clipping bug, and pass **`--no_poly_g`** to preserve the pipeline's earlier default. Updating the version alone does not preserve that default.
- Give Chromap technical runs from the same biological library a shared **`RGLB`**, while retaining distinct **`RGID`** values. Picard can then mark duplicates across those runs. This is independent of trimming.
- Preserve PRESEQ's original command options on retry and append **`-defects`**. The updated wrapper calls `String.join` with no elements, producing an empty string and dropping `-bam`, `-verbose` and `-seed 1`. Record the local module patch in `modules.json`.

- Update GFFREAD to **0.12.9**, using the packaged Bioconda container and matching Conda declaration. This fixes an upstream annotation-conversion regression affecting valid GFF scores. The pipeline GFF-to-GTF command is preserved, and the local module patch is registered. The historical HIV runs use GTF input and do not exercise this branch.

The PRESEQ correction preserves the module's new retry policy. The old ignore-to-new-retry policy, which can disable the estimate's defect check, is a separate documented behavior change; successful completion does not establish the validity of every estimate.

## Validation

- Four full-input HIV workflows succeeded: BWA and Chromap on reviewed dev `73155d` and the frozen two-fix proposal `a9b588`. Each processed **87,089,251 pairs** across four libraries.
- The fixed BWA library-consensus matrix exactly matches the archived original dev matrix: **74,895 regions and all 299,580 count cells**. The unfixed matrix matches the archived PR448 matrix. This is an archived-matrix comparison, not a fresh run of the immediate pre-448 parent; historical FASTQ byte identity is unestablished.
- With identical real alignments from two resequenced yeast libraries, pre-448 and fixed duplicate marking agrees exactly. Unfixed dev retains **2,553 and 1,865 additional duplicate pairs**.
- Three yeast core workflow runs restore individual peak coordinates and all **1,934 consensus regions**. One count cell differs by one, traced to an alignment/duplicate-representative change. Optional QC exclusions were applied equally. The HIV cohort has one technical run per library and cannot isolate the library-ID defect.
- On one full cleaned human HIV library BAM (285,036,036 bytes), MACS3 3.0.1 and 3.0.4 produce exactly the same **30,438 broad peaks**, every broadPeak field and every non-comment XLS row, using the archived default broad/BAMPE/no-control arguments. The 3.0.4 replay also matches archived cloud outputs. Only gappedPeak display fields 7–8 differ between versions. Four offline comparison tests pass; this covers the emitted precision and one library/configuration, not every mode or biological equivalence.
- `python3 tests/regressions/preseq_retry.py` exercises the actual PRESEQ wrapper: normal success, retry recovery and an intentionally broken mutation. Its executable is an argument-checking stand-in, so it does not validate numerical estimates.
- `python3 tests/regressions/module_patch_roundtrip.py`, in an nf-core/tools environment, reapplies the registered patches to the exact `modules.json` revisions and compares the installed files byte for byte. The previous GFFREAD patch fails this check; the regenerated patch passes with nf-core/tools 4.0.3 and 4.1.0, alongside the PRESEQ and Trim Galore controls. Snapshot JSON values are unchanged; its final newline matches the upstream format.
- Changed-file pre-commit checks, including Nextflow lint, pass. Full nf-core/tools 4.1.0 lint reports 13 failures: five template differences and eight container-config differences. A separate run on immutable `a9b588` reproduces the same 13 failures. These remain documented baseline debt; full lint is not reported as passing. Its formatter setup also encountered a local npm/prek error, while the separately executed repository formatter hook passes.

All nine existing GFFREAD module tests pass locally in `--ci` mode without snapshot updates, using nf-test 0.9.5, Nextflow 26.04.6 and the upstream module helper nft-utils 0.0.9. Reviewed snapshot changes are version values, GFF3 comment headers and run metadata; GTF, FASTA, BED and stub hashes are unchanged. The pipeline's separate test configuration still pins nft-utils 0.0.3. Docker amd64 was exercised; matching Conda declarations and SIF availability were checked, while Conda resolution and SIF execution were not.

The historical full-data runs remain pinned to `a9b588`; they have not been rerun or relabeled as testing the later PRESEQ or GFFREAD corrections. [Reproducibility repository](https://github.com/NicoDeVeaux/atacseq-modules-update-verification). The full output-tree archive is being prepared in an unpublished draft release.

## Broader review and opening gates

The [module audit](docs/module-update-audit.md) accounts for all **267 changed module/subworkflow files**, **70 components**, **16 changed declared dependency versions**, and subsequent histories for **36 shared module paths and seven shared subworkflows**.

- [x] Reproduce the two original regressions and validate their mitigations.
- [x] Scan every declared version change and later merged module fix.
- [x] Reproduce PRESEQ retry option loss and exercise its correction.
- [x] Classify the tested SAMtools and Picard candidates as inherited defects using old/selected binary controls.
- [x] Classify the tested deepTools boundary behavior as inherited.
- [x] Identify the GFFREAD conversion regression and prepare its 0.12.9 mitigation.
- [x] Validate the packaged GFFREAD conversion, all nine existing module tests and the reviewed snapshots.
- [x] Finish the bounded alignment/peak review: real mappings, QC partitions, SAM metadata and display fields are compared separately. The human control below covers one broad/no-control library; SE/ALT/repetitive-data equivalence remains unestablished.
- [ ] Finish and link the downloadable computation/output archive.
- [ ] Record final branch validation and update this description before opening the nf-core PR.

Chromap 0.3.2 includes a mapping-pair correction. Local tests on the real residual pair reproduce an alignment difference independently of trimming and read-group assignment. Complete alignment identity with 0.2.4 is therefore not promised. Later upstream fixes are leads to investigate, not proof that every affected version causes a new pipeline regression. This branch does not blanket-upgrade other tools or resolve the separate mixed SE/PE and MultiQC issues.

This is written by GPT-6 (Codex) on Nico DeVeaux's behalf.
