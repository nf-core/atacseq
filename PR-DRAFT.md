# Draft PR: fix module-update regressions after #448

Target: `nf-core/atacseq:dev`  
Related issue: [#473](https://github.com/nf-core/atacseq/issues/473)  
Status: versioned proposal on the fork; module review, patch-lifecycle repair and archive preservation complete. Final targeted pipeline checks remain pending. See the [qualification plan](docs/pr473-qualification.md).

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
- A separate local execution uses the actual PRESEQ 3.2.0 package, current module, production configuration selectors and two existing fixed-arm yeast merged, duplicate-marked BAMs. Both tasks succeed on attempt one with `-verbose -bam -seed 1 -pe`, producing 10,000 finite curve rows each. Input hashes remain unchanged and match the archive. No real retry occurs; this establishes execution on these inputs, not estimation accuracy or reliability of extreme extrapolation.
- `python3 tests/regressions/module_patch_roundtrip.py`, in an nf-core/tools environment, checks registered application, canonical regeneration, and forward/reverse replay against the exact `modules.json` revisions. The earlier forward-only repair did not survive regeneration. The extended check fails before the new correction and passes with nf-core/tools 4.0.3 and 4.1.0, alongside PRESEQ and Trim Galore controls. Fresh Nextflow 25.10.4 snapshots preserve all eight scientific content values and the upstream EOF format.
- Nextflow 26.04.6 lint reports zero errors and three existing warnings across 131 files. Full nf-core/tools 4.0.3 lint passes in the official Linux/amd64 container on both repair `81bea6` and upstream dev `73155d`: 515 passes, four ignored, 127 warnings and zero failures each. The corrected GFFREAD module separately has 47 passing lint checks, one newer-module warning and no failures. The four earlier macOS template mismatches were reproduced on both source trees and resolved by normal Linux formatting/rendering without source changes or weakened checks. The container used Python 3.14.4 and Nextflow 25.10.4; its unmodified hook resolver selected Node 27 alpha/npm 12.1, so this is local qualification rather than exact historical CI dependency parity.

All nine existing GFFREAD module tests pass locally in `--ci` mode without snapshot updates on Nextflow 25.10.4 and, in the earlier run, 26.04.6, using nf-test 0.9.5 and upstream helper nft-utils 0.0.9. The new snapshot metadata comes from the actual 25.10.4 run. Reviewed scientific contents are unchanged. The pipeline's separate test configuration still pins nft-utils 0.0.3 and ignores vendored module tests. Docker amd64 was exercised; matching Conda declarations and SIF availability were checked, while Conda resolution and SIF execution were not.

The historical full-data runs remain pinned to `a9b588`; they have not been rerun or relabeled as testing the later PRESEQ or GFFREAD corrections. [Public reproducibility fixtures](https://github.com/NicoDeVeaux/atacseq-modules-update-verification) remain available. Validation data are preserved in verified private Seqera S3 archives. Local copies are retained until #473 is resolved.

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
- [x] Finish verified S3 archive closeout and retain the public reproducibility links.
- [ ] Run the two small final-SHA GTF/GFF3 pipeline checks with real PRESEQ described in the qualification plan.
- [ ] Record final branch validation and update this description before opening the nf-core PR.

Chromap 0.3.2 includes a mapping-pair correction. Local tests on the real residual pair reproduce an alignment difference independently of trimming and read-group assignment. Complete alignment identity with 0.2.4 is therefore not promised. Later upstream fixes are leads to investigate, not proof that every affected version causes a new pipeline regression. This branch does not blanket-upgrade other tools or resolve the separate mixed SE/PE and MultiQC issues.

Earlier validation narrative attribution: GPT-6 (Codex). This review update is written by an AI agent (model name unavailable) on Nico DeVeaux's behalf.
