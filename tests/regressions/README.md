# Focused module regression checks

Run `python3 tests/regressions/preseq_retry.py` with Python 3 and Nextflow installed. The check uses offline mode, so its Nextflow distribution must already be cached. `--nextflow-jar /path/to/nextflow-one.jar` selects an explicit distribution; `--output /empty/directory` retains task commands, logs, traces and results in a chosen directory.

The check executes the actual PRESEQ module on the local executor without containers or Conda. Its stand-in executable records arguments, simulates an initial complexity-model failure, and succeeds only when original options and the retry flag are present. It checks normal success, retry recovery and an intentionally broken mutation. It does not validate scientific estimates or the complete ATAC pipeline.

The real-data comparisons and wider review are recorded in `PR-DRAFT.md` and `docs/module-update-audit.md`.

Run `python3 tests/regressions/gffread_conversion.py` for the actual packaged GFFREAD conversion check. It requires Docker and tests feature coordinates, identifiers and the observed fixed runtime on a benign annotation. Use `--docker-context colima` on the local Mac. It does not reproduce historical failures or prove all annotation inputs equivalent.

Run `python3 tests/regressions/module_patch_roundtrip.py` in a Python environment containing nf-core/tools (validated with the repository's 4.0.3 version and with 4.1.0). It fetches the exact upstream revisions in `modules.json`, applies the registered GFFREAD and PRESEQ patches, regenerates them with nf-core's default diff settings, and checks forward and reverse replay byte for byte. Trim Galore is an unpatched control. `--cache /path/to/cache` reuses fetched files and `--output /empty/directory` retains results. This is a manual regression check, separate from the pipeline's nf-test selection.
