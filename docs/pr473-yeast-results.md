# Matched native yeast qualification

The two fresh workflows at source `4def83` pass the targeted execution and numerical gates on Nextflow 26.04.6. This covers the final branch's native GTF/GFF3 pipeline paths with real PRESEQ, alongside the earlier module and causal regression controls.

| Check                      | Observed result                                                                                                                                     |
| -------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| Execution                  | GTF: 92 completed tasks, about 27 minutes. GFF3: 93 completed tasks, about 22 minutes. All attempts were one; neither run used resume.              |
| Read groups and containers | Correct four technical read groups and two physical-library groups; actual image digests match for common processes.                                |
| Marked/cleaned BAMs        | Exact read and retained-fragment multisets agree. D28 retains 282,467 fragments; D28P retains 283,754.                                              |
| Peaks and counts           | Exact coordinates: 1,635 D28 peaks, 1,431 D28P peaks and 1,934 consensus intervals. All 3,868 keyed count cells agree.                              |
| PRESEQ                     | Both libraries succeed in both arms with original arguments and 10,000 finite rows each. Matched curves are byte-identical. No real retry occurred. |
| Annotation path            | GUNZIP_GFF and GFFREAD 0.12.9 run once in GFF3; GTF runs no GFFREAD task. Recognized RNA geometry and normalized IDs pass.                          |

These values were measured on the new native Ensembl runs, not imported from the older UCSC oracle. The checker has 34 benign controls covering missing/stale evidence, scientific-parameter differences, actual image-digest changes, MultiQC and incomplete execution. Command, version, attempt and digest evidence is bound to the actual task records.

The first GFF startup failed before scientific tasks: the schema accepts `.gff`/`.gff.gz`, while the official Ensembl filename ends in `.gff3.gz`. The successful run uses a byte-identical `.gff.gz` alias. Validation stayed enabled, application source and annotation content stayed unchanged, and the rejected startup is retained in the report. Replays must use the supported suffix.

Default GFFREAD arguments omit exactly 91 non-RNA transposable-element transcript features: the converted annotation has 7,036 transcripts/7,416 exons, versus 7,127/7,507 in native GTF. Namespace differences and the omission are classified separately from numerical core equivalence; the branch does not claim complete annotation byte identity.

This is targeted software qualification on two genuine resequenced yeast libraries, not a differential-biology study, proof of PRESEQ estimation accuracy or all-mode certification. Conda/SIF execution and the deliberately excluded optional QC modes remain outside this campaign. The existing mutation test remains the retry witness.

[Public scientific outputs, manifests, offline verifier and report](https://github.com/NicoDeVeaux/atacseq-modules-update-verification/tree/da54c0ee85bcb11ce088e1d974550799fe13edfb/verification/reviews/2026-11-02-pr473-yeast). Every published blob and the joined compressed archive received anonymous SHA256 readback. The portable checker also passes on the genuine data extracted from the sealed archive. This closeout changes review documentation only after tested source `4def83`; runtime, tests, configuration and container declarations are unchanged. Final Nextflow lint again reports 131 files, zero errors and three existing warnings. Historical full-depth evidence keeps its original source pins. All local validation copies remain retained until #473 is resolved. No production change or upstream PR was made.

This is written by an AI agent (model name unavailable) on Nico DeVeaux's behalf.
