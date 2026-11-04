# CI regression checks

`preseq_retry.py` runs in the blocking **PRESEQ retry regression** job in `.github/workflows/nf-test.yml`. The existing `confirm-pass` gate depends on it, independently of changed-test selection.

The check runs the actual PRESEQ Nextflow module with a stand-in executable: normal success, failure followed by successful retry, and a negative control whose configured arguments are deliberately lost after the first attempt. It verifies that verbosity, BAM mode and the fixed seed survive retry, and that `-defects` is added only on retry. It does not assess PRESEQ estimate quality.

For local reproduction, cache Nextflow with `nextflow -version`, then run `python3 tests/regressions/preseq_retry.py`. `--nextflow-jar /path/to/nextflow-one.jar` selects a cached distribution; `--output /empty/directory` retains commands, logs, traces and results. CI uploads that fixture directory after success or failure.

The opt-in GFFREAD conversion and module-patch replay checks live in the [supporting-data repository](https://github.com/NicoDeVeaux/atacseq-modules-update-verification/tree/bca333aa451587a3c649224b0edbfcabc618f647/verification/checks/issue-473). They require an explicit `--repo` pointing at the pipeline checkout and are not pipeline CI tests.

The blocking **Trimming and duplicate-library regressions** job also runs:

- `trimgalore_poly_g.py`: the production Trim Galore module and configuration must preserve all 32 synthetic read pairs, including legitimate R1 poly-G and R2 5′ poly-C. Checks compare read IDs, sequences and qualities directly.
- `duplicate_libraries.py`: the production Picard AddOrReplaceReadGroups and MarkDuplicates modules must mark one resequenced fragment as duplicate while retaining a coincident fragment from another biological library. All five pairs must survive. This isolates library labels and marking; it does not exercise full Chromap alignment or workflow grouping.

Both require Nextflow 25.10.4 and cached Docker images, listed in the CI workflow. Run each script with `--output /empty/directory`; `--repo /pipeline/checkout` selects another source version and `--docker-context NAME` selects Docker. `NXF_VER=25.10.4` selects Nextflow. CI retains source hashes, commands, traces and results for seven days. Unit controls run with `python3 -B -m unittest discover --start-directory tests/regressions --pattern 'test_*.py'`.
