# Change for review — fix `ch_fasta_fai` reference-channel broadcast

_Author: Co-Scientist · Date: 2026-07-05 · Branch: `nf-core-modules-update` (fork `NicoDeVeaux/atacseq`)_
_Reviewer: please focus on the single change described below and the "How to verify" section._

---

## TL;DR

One-line fix in `workflows/atacseq.nf`: append `.first()` to the `ch_fasta_fai`
channel so the genome reference is a **value channel** and broadcasts to every
sample. Without it, the pipeline silently processes **only the first sample**.

```diff
     ch_fasta_fai = ch_fasta
         .combine(ch_fai)
         .map { fasta, fai -> [ [:], fasta, fai ] }
+        .first()
```

(+ an explanatory comment above the block.)

This is a **correctness** fix, not cosmetic. It has no effect on single-sample
runs and is required for correct multi-sample runs (i.e. all real ATAC-seq
datasets and the `test_full` profile).

---

## Where the bug came from

This branch combines two efforts: a strict-syntax migration and an
`nf-core modules update`. The module update changed the input contract of the
shared alignment/sort/markduplicates components: they now take a **single
combined reference tuple** `[ meta, fasta, fai ]` instead of separate `fasta`
and `fai` inputs. To feed them, new glue code was added in `workflows/atacseq.nf`
(~line 169):

```nextflow
ch_fasta_fai = ch_fasta
    .combine(ch_fai)
    .map { fasta, fai -> [ [:], fasta, fai ] }
```

`ch_fai` is derived from `SAMTOOLS_FAIDX.out.fai` — a **process output**, which
in DSL2 is a **queue** channel (here, exactly one element). `combine` on a queue
returns a queue, and `.map` preserves that. So `ch_fasta_fai` was a
**single-element queue channel**.

### Why that breaks multi-sample runs

When a Nextflow process receives one multi-item queue input (the per-sample
BAM/reads) and one single-item queue input (the reference), it executes only
**`min(N, 1) = 1`** time — the run stops as soon as the shorter queue is
exhausted. Value channels, by contrast, are broadcast (re-consumed) for every
item of the driving queue. This is the canonical Nextflow reference-channel
footgun, and the standard remedy is to make the reference a value channel via
`.first()` (or `.collect()`).

`.first()` converts the one-element queue into a value channel, restoring the
broadcast so all N samples are processed.

### Why `dev` did not have this bug

On `dev`, the reference was passed as `ch_fasta` — built once as
`channel.value(file(params.fasta))`, i.e. already a **value** channel — directly
into the old subworkflows. The value-channel property was never lost. The module
update's new `combine`-based tuple wiring dropped that property, and the fix
restores it.

---

## Blast radius (what was affected before the fix)

`ch_fasta_fai` is consumed as the reference input by, among others:

- `FASTQ_ALIGN_BWA` / `FASTQ_ALIGN_BOWTIE2` / `FASTQ_ALIGN_CHROMAP` → `BWA_MEM` etc.
- `BAM_SORT_STATS_SAMTOOLS` → `SAMTOOLS_SORT`, `BAM_STATS_SAMTOOLS`
- `BAM_MARKDUPLICATES_PICARD` → `PICARD_MARKDUPLICATES`
- `subworkflows/local/bam_filter_bamtools.nf`, `subworkflows/local/align_star.nf`

None of these re-broadcast the reference internally (they only `.map` it), so the
single-element-queue property propagated everywhere. Net effect before the fix:
**only the first sample would flow through alignment → sort → dedup → downstream**;
remaining samples were silently dropped. The run could still report `SUCCEEDED`
with biologically wrong/incomplete output, or fail later at a multi-sample
grouping/consensus step.

---

## Evidence (independently reproducible)

Minimal reproduction of exactly this wiring (3 samples + a `combine`-derived
one-element reference), counting process executions:

| Wiring                                     | Executions (3 samples in) |
| ------------------------------------------ | ------------------------- |
| `combine(...).map{...}` (pre-fix)          | **1** ❌                  |
| `combine(...).map{...}.first()` (post-fix) | **3** ✅                  |

Reproduce with a ~15-line script: a process taking `val sample` + `val ref`,
driven by `channel.of('s1','s2','s3')` and
`channel.value('FASTA').combine(channel.of('FAI')).map{...}` with and without
`.first()`; count how many times it runs.

---

## Verification

- `nextflow lint .` → **0 errors** (6 pre-existing single-emit style warnings,
  unrelated to this change).
- **Still pending (recommended before merge):** a multi-sample `test_full` run on
  the fixed commit must reach `SUCCEEDED` and show alignment/sort/markduplicates
  tasks equal to the sample count (not 1). The in-flight run
  `intergalactic_poitras` (`1kEThYNVCX9ZOK`) was launched on the _pre-fix_ commit
  and is therefore invalid — cancel it and relaunch on the commit containing this
  fix.

---

## Related item NOT changed (flagged for reviewer judgement)

`ch_fasta` itself is a value channel for a plain FASTA but a **one-element queue**
when the input genome is gzipped (`GUNZIP_FASTA(...).gunzip.map{...}` in
`subworkflows/local/prepare_genome.nf`). `subworkflows/local/bam_shift_reads.nf`
consumes `ch_fasta` directly (not `ch_fasta_fai`) and would hit the same
single-sample cap **only when the reference is gzipped**. `test_full` uses an
unzipped iGenomes FASTA, so it does not trigger there, and I left it unchanged to
keep this PR minimal. A sturdier belt-and-suspenders alternative would be to make
`ch_fasta`/`ch_fai` value channels at the `PREPARE_GENOME` output boundary (e.g.
`.first()`), which would also make the `ch_fasta_fai` fix redundant. Reviewer's
call on whether to fold that in.

---

## What I explicitly verified as correct (no change needed)

- **Topic-channel version handling** — 53 modules emit
  `tuple val("${task.process}"), val('tool'), eval(...), topic: versions`; zero
  still emit `versions.yml`; the collation
  `channel.topic('versions').map { process, tool, version -> ... }.groupTuple().collectFile(...)`
  matches the emitted tuple shape.
- **Module-signature reconciliations** — `SAMTOOLS_SORT` (3 args incl.
  `index_format` at all call sites), `GFFREAD` (2 args, `[]` for optional fasta),
  `PICARD_MARKDUPLICATES` (2 args, fasta+fai combined) all match their updated
  module input blocks.
- **Strict-syntax rewrites** — `Channel.` → `channel.`, explicit closure params,
  thinned `NFCORE_ATACSEQ` entry in `main.nf`: no issues.
