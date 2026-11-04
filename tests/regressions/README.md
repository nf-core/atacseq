# SAF boundary regression

The blocking **SAF coordinate boundary regression** job runs `saf_boundary.py` against the converter extracted from `modules/local/macs3_consensus.nf`, using real Subread 2.0.1 and Samtools runtimes pinned by digest. It verifies zero and nonzero BED starts, feature identity, and assignment of three preserved single-base alignments. The read immediately before the corrected interval must be excluded; two positive controls must remain counted.

Cache the two images listed in `.github/workflows/nf-test.yml`, then run `python3 -B tests/regressions/saf_boundary.py --output /empty/directory`. Use `--repo /pipeline/checkout` to select another source version and `--docker-context NAME` to select Docker. Run checker controls with `python3 -B -m unittest discover --start-directory tests/regressions --pattern 'test_saf_boundary.py'`.

CI retains source hashes, runtime identities, commands and results for seven days. This test executes the extracted converter and real counting, not the complete MACS3 consensus workflow. Existing nf-test assertions separately compare BED/SAF pairs from full default-profile runs.
