#!/bin/bash
# The filesystem roots below were specific to the machines this pipeline
# was run on. They are replaced by /path/to/machine and /path/to/home
# placeholders for deposition; repoint them at your own storage before
# running this script.
# Seed the fan-out with the sample_08 output already produced by smoke_a, so
# array index 8 is idempotent-skipped rather than recomputed.  A hard link, not
# a copy: same filesystem, same inode, so the bytes are provably identical and
# no extra 3.9 GB is spent.
set -euo pipefail
ROOT=/path/to/home/graphst_repro
SRC=$ROOT/run/smoke_a/results/sample_08_graphst_output.h5ad
DST=$ROOT/run/full/results/sample_08_graphst_output.h5ad
mkdir -p "$(dirname "$DST")"
[ -s "$SRC" ] || { echo "FATAL: smoke_a output missing"; exit 1; }
[ -e "$DST" ] && { echo "already present: $DST"; exit 0; }
ln "$SRC" "$DST"
echo "linked  $DST"
echo "inode   $(stat -c %i "$SRC") == $(stat -c %i "$DST")"
echo "bytes   $(stat -c %s "$DST")"
