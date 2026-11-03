# CI regression check

`preseq_retry.py` runs in the blocking **PRESEQ retry regression** job in `.github/workflows/nf-test.yml`. The existing `confirm-pass` gate depends on it, independently of changed-test selection.

The check runs the actual PRESEQ Nextflow module with a stand-in executable: normal success, failure followed by successful retry, and a negative control whose configured arguments are deliberately lost after the first attempt. It verifies that verbosity, BAM mode and the fixed seed survive retry, and that `-defects` is added only on retry. It does not assess PRESEQ estimate quality.

For local reproduction, cache Nextflow with `nextflow -version`, then run `python3 tests/regressions/preseq_retry.py`. `--nextflow-jar /path/to/nextflow-one.jar` selects a cached distribution; `--output /empty/directory` retains commands, logs, traces and results. CI uploads that fixture directory after success or failure.

The opt-in GFFREAD conversion and module-patch replay checks live in the [supporting-data repository](https://github.com/NicoDeVeaux/atacseq-modules-update-verification/tree/bca333aa451587a3c649224b0edbfcabc618f647/verification/checks/issue-473). They require an explicit `--repo` pointing at the pipeline checkout and are not pipeline CI tests.
