# Qualification before the upstream #473 PR

The original trimming and Chromap library-ID acceptance checks are supported by the frozen `a9b588` runs. Later runtime changes affect GFFREAD and PRESEQ, which those GTF-input, PRESEQ-disabled workflows did not execute. Keep that distinction when describing the final branch.

The 2026-11-02 review found and repaired one maintenance issue: applying the GFFREAD patch worked, but canonical nf-core regeneration produced a patch that could not replay. The regression now checks registered application, default regeneration, and forward/reverse replay. It fails before the repair and passes with nf-core/tools 4.0.3 and 4.1.0. Fresh benign snapshots were recorded with Nextflow 25.10.4; all eight scientific snapshot contents are unchanged. The strict module suite passes all nine tests without updating snapshots.

Full nf-core/tools 4.0.3 lint passes in the official Linux/amd64 container on repair `81bea6` and upstream dev `73155d`: 515 passes, four ignored, 127 warnings and zero failures each. Linux regenerated all three committed logos byte for byte and ran formatting normally; the four local macOS template mismatches were environment/toolchain artifacts. No checks or security settings were weakened. Nextflow 26.04.6 lint reports zero errors and three existing single-emit warnings across 131 files. The corrected GFFREAD module passes 47 lint checks with one newer-module warning. The Linux environment used Python 3.14.4/Nextflow 25.10.4 and its ordinary hook resolver selected Node 27 alpha/npm 12.1; this is local qualification, not exact historical CI dependency parity.

A separate local check executes the actual PRESEQ 3.2.0 module with production configuration selectors on both fixed-arm marked yeast BAMs. Both tasks succeed on attempt one, retain the production arguments, and emit 10,000 finite curve rows each. Their input hashes remain unchanged. No real retry occurs; this adds packaged-tool execution evidence without validating estimation accuracy or replacing a full pipeline check.

## Executed workflow checks and reproduction contract

The two matched workflows completed at immutable source `4def83`, with Chromap and real PRESEQ enabled. Reproduce them using the following contract. Use the genuine D28/D28P technical resequencing inputs from PRJNA844120:

- D28: SRR19479022 and SRR19479023, two runs from the same physical library.
- D28P: SRR19479020 and SRR19479021, two runs from a different physical library.
- Start with the existing first 250,000 complete pairs per run, preserving their recorded byte hashes and samplesheet grouping.

| Arm  | Annotation                       | Required path                                                                 |
| ---- | -------------------------------- | ----------------------------------------------------------------------------- |
| GTF  | Ensembl R64-1-1 release-115 GTF  | No GFFREAD task; real PRESEQ results for both libraries                       |
| GFF3 | Matching release-115 native GFF3 | GUNZIP_GFF and GFFREAD 0.12.9 execute; real PRESEQ results for both libraries |

Use the matching Ensembl FASTA in both arms. Stage the official GFF3 bytes under a `.gff.gz` alias because the current schema rejects the `.gff3.gz` filename suffix. Explicitly ignore iGenomes, leave `genome` unset, set `mito_name=Mito`, and ensure the GFF arm has effective `gtf=null`. Freeze all other scientific settings and container identities before launch. Validate annotation contigs and interval bounds against the complete FASTA first.

Require four distinct technical read groups, one shared library ID per physical library, two separate biological-library groups, and both samples in downstream results. Compare retained fragments, peak coordinates, consensus intervals and coordinate/sample-keyed count cells exactly between the two arms. Validate converted annotation identifiers and coordinates separately; native GTF/GFF3 representation differences are not automatically regressions. These references do not inherit the old UCSC 1,934-region oracle.

Inspect actual PRESEQ commands, versions, attempts and parseable estimates. Original `-verbose -bam -seed 1 -pe` options must survive any retry; `-defects` is added only on retry. A normal run need not retry. The existing mutation check establishes retry argument handling; successful completion alone does not establish numerical model validity. Record insufficient-data or model failures rather than suppressing them.

Prior bounded yeast runs took about 35–37 minutes per arm; allow roughly 30–60 minutes per arm plus staging, queues and cold images. This is an estimate, not a dollar quote or a measured fresh cloud run. Another four-arm, 87-million-pair campaign is unnecessary for these unchanged core paths. Full human reruns remain an optional stronger final-SHA qualification.

## Review evidence and opening gates

Reuse successful nf-core/nf-test CI when it tested the same SHA, or when every relevant source, test, configuration and container input is demonstrably unchanged. Record that correspondence and confirm the relevant suites actually ran. An all-skipped aggregate success alone is insufficient. Always run Nextflow lint. Keep module/subworkflow checks separate from pipeline execution and numerical output checks.

Keep upstream PRESEQ/GFFREAD module work visible. Registered local patches are deliberate temporary integration points; document their pins and lifecycle checks rather than implying they are unmodified upstream modules. Do not broaden this PR to unrelated template, MultiQC, mixed SE/PE or inherited aligner behavior.

The two real-sample arms have completed and passed the targeted gates; see the [measured results and limits](pr473-yeast-results.md). Preserve all local validation data until [nf-core/atacseq#473](https://github.com/nf-core/atacseq/issues/473) is resolved. Historical full-data storage is private Seqera S3. The new yeast campaign publishes scientific outputs and its report to the public data repository linked in the measured results. The upstream PR has intentionally not been opened.

References: [Evangelos's release-planning guidance](https://github.com/vagkaratzas/skills/blob/4f3aebbc0dc776c93b266cef51f74f92aaeac5fa/skills/plan-nf-core-releases/SKILL.md), [PRJNA844120](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA844120), [Ensembl FASTA](https://ftp.ensembl.org/pub/release-115/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz), [GTF](https://ftp.ensembl.org/pub/release-115/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.115.gtf.gz), [GFF3](https://ftp.ensembl.org/pub/release-115/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.115.gff3.gz).

This is written by an AI agent (model name unavailable) on Nico DeVeaux's behalf.
