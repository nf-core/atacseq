#!/usr/bin/env bash
#
# Unit tests for bin/featurecounts_merge.sh — pure bash, no Docker, milliseconds.
#
# The consensus-peak quantification runs featureCounts once per library type
# (single-end / paired-end) and merges the resulting count tables. These tests
# pin the merge logic directly: exact column-bind of a mixed SE+PE pair, and the
# single-file pass-through used when a cohort is all-SE or all-PE.
#
# Run: tests/bin/test_featurecounts_merge.sh
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
merge="$here/../../bin/featurecounts_merge.sh"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

fail=0
pass=0

# assert_eq NAME EXPECTED_FILE ACTUAL_FILE
assert_eq() {
    local name=$1 exp=$2 act=$3
    if diff -u "$exp" "$act" >/dev/null; then
        echo "ok   - $name"
        pass=$((pass + 1))
    else
        echo "FAIL - $name"
        echo "----- diff (expected vs actual) -----"
        diff -u "$exp" "$act" || true
        echo "-------------------------------------"
        fail=$((fail + 1))
    fi
}

# --- fixtures --------------------------------------------------------------
# Both tables share an identical annotation block (same SAF); they differ only
# in their sample count columns. The PE table intentionally lists its data rows
# in a DIFFERENT order to prove the merge keys on Geneid, not on row position.

printf '%s\n' \
'# Program:featureCounts v2.1.1; Command:"featureCounts" "-F" "SAF"' \
$'Geneid\tChr\tStart\tEnd\tStrand\tLength\tT100_SE.bam\tT150_SE.bam' \
$'peak_1\tI\t1\t100\t+\t100\t11\t12' \
$'peak_2\tII\t5\t205\t+\t201\t21\t22' \
$'peak_3\tIII\t9\t309\t+\t301\t31\t32' \
> "$tmp/se.featureCounts.tsv"

printf '%s\n' \
'# Program:featureCounts v2.1.1; Command:"featureCounts" "-F" "SAF" "-p"' \
$'Geneid\tChr\tStart\tEnd\tStrand\tLength\tT0_PE.bam\tT15_PE.bam' \
$'peak_3\tIII\t9\t309\t+\t301\t131\t132' \
$'peak_1\tI\t1\t100\t+\t100\t111\t112' \
$'peak_2\tII\t5\t205\t+\t201\t121\t122' \
> "$tmp/pe.featureCounts.tsv"

# === Case 1: mixed SE + PE merge (SE table first) ==========================
# Annotation and row order come from the first (SE) file; sample columns are
# appended SE-then-PE; PE counts are matched to rows by Geneid.
printf '%s\n' \
'# Program:featureCounts (merged single-end and paired-end libraries)' \
$'Geneid\tChr\tStart\tEnd\tStrand\tLength\tT100_SE.bam\tT150_SE.bam\tT0_PE.bam\tT15_PE.bam' \
$'peak_1\tI\t1\t100\t+\t100\t11\t12\t111\t112' \
$'peak_2\tII\t5\t205\t+\t201\t21\t22\t121\t122' \
$'peak_3\tIII\t9\t309\t+\t301\t31\t32\t131\t132' \
> "$tmp/expected_mixed.tsv"

"$merge" "$tmp/out_mixed.tsv" "$tmp/se.featureCounts.tsv" "$tmp/pe.featureCounts.tsv"
assert_eq "mixed SE+PE merge column-binds on Geneid" "$tmp/expected_mixed.tsv" "$tmp/out_mixed.tsv"

# === Case 2: all-SE cohort -> single-file pass-through =====================
# Empty PE branch means featureCounts runs once; the merge must pass the single
# table through unchanged (except the normalised comment line), order preserved.
printf '%s\n' \
'# Program:featureCounts (merged single-end and paired-end libraries)' \
$'Geneid\tChr\tStart\tEnd\tStrand\tLength\tT100_SE.bam\tT150_SE.bam' \
$'peak_1\tI\t1\t100\t+\t100\t11\t12' \
$'peak_2\tII\t5\t205\t+\t201\t21\t22' \
$'peak_3\tIII\t9\t309\t+\t301\t31\t32' \
> "$tmp/expected_se_only.tsv"

"$merge" "$tmp/out_se_only.tsv" "$tmp/se.featureCounts.tsv"
assert_eq "all-SE single-file pass-through" "$tmp/expected_se_only.tsv" "$tmp/out_se_only.tsv"

# === Case 3: all-PE cohort -> single-file pass-through =====================
# Same pass-through, preserving the PE file's own (unsorted) row order.
printf '%s\n' \
'# Program:featureCounts (merged single-end and paired-end libraries)' \
$'Geneid\tChr\tStart\tEnd\tStrand\tLength\tT0_PE.bam\tT15_PE.bam' \
$'peak_3\tIII\t9\t309\t+\t301\t131\t132' \
$'peak_1\tI\t1\t100\t+\t100\t111\t112' \
$'peak_2\tII\t5\t205\t+\t201\t121\t122' \
> "$tmp/expected_pe_only.tsv"

"$merge" "$tmp/out_pe_only.tsv" "$tmp/pe.featureCounts.tsv"
assert_eq "all-PE single-file pass-through" "$tmp/expected_pe_only.tsv" "$tmp/out_pe_only.tsv"

# === Case 4: usage error on too few arguments =============================
if "$merge" "$tmp/out_none.tsv" >/dev/null 2>&1; then
    echo "FAIL - merge should exit non-zero with no input files"
    fail=$((fail + 1))
else
    echo "ok   - usage error when no input tables given"
    pass=$((pass + 1))
fi

# --- summary ---------------------------------------------------------------
echo
echo "featurecounts_merge.sh: $pass passed, $fail failed"
[ "$fail" -eq 0 ]
