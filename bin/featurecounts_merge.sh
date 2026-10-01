#!/usr/bin/env bash
# Merge multiple featureCounts count tables into a single matrix.
#
# The consensus-peak quantification runs featureCounts once per library type
# (single-end / paired-end), so each input table shares an identical annotation
# block (Geneid, Chr, Start, End, Strand, Length) computed from the same SAF, and
# differs only in its per-sample count columns. This column-binds those sample
# columns back together, keyed on Geneid, and reproduces the featureCounts output
# layout expected downstream by deseq2_qc.r:
#
#   line 1        : a "# Program:featureCounts" comment (skipped via read.delim skip=1)
#   line 2        : header  Geneid  Chr  Start  End  Strand  Length  <sample cols...>
#   remaining rows: annotation columns 1-6 followed by one count per sample
#
# With a single input table this is an order-preserving pass-through.
#
# Usage: featurecounts_merge.sh OUTFILE INPUT1 [INPUT2 ...]
set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "Usage: $(basename "$0") OUTFILE INPUT1 [INPUT2 ...]" >&2
    exit 1
fi

out=$1
shift

awk '
    BEGIN { FS = OFS = "\t" }

    # Skip the leading "# Program:featureCounts ..." comment of every input file.
    FNR == 1 { fidx++; next }

    # Header row: keep the 6 annotation columns from the first file only,
    # then append every input file`s sample columns (7..NF) in argument order.
    FNR == 2 {
        if (fidx == 1) { hdr = $1; for (i = 2; i <= 6; i++) hdr = hdr OFS $i }
        for (i = 7; i <= NF; i++) hdr = hdr OFS $i
        next
    }

    # Data rows: index on Geneid (column 1). Preserve the first file`s row order
    # and annotation; append sample counts from each file for the matching Geneid.
    {
        key = $1
        if (fidx == 1) {
            order[++n] = key
            a = $1; for (i = 2; i <= 6; i++) a = a OFS $i; ann[key] = a
            v = "";  for (i = 7; i <= NF; i++) v = v OFS $i; val[key] = v
        } else {
            for (i = 7; i <= NF; i++) val[key] = val[key] OFS $i
        }
    }

    END {
        print "# Program:featureCounts (merged single-end and paired-end libraries)"
        print hdr
        for (j = 1; j <= n; j++) print ann[order[j]] val[order[j]]
    }
' "$@" > "$out"
