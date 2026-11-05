#!/bin/bash
# Build outputs/gsea_13types: the Hallmark tables of twelve cell types
# (paths.SOUND_GSEA_DIR) plus the neutrophil tables written by run.sh.
# Run after run.sh: bash merge_thirteen.sh
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
MOD="$(dirname "$HERE")"

# The twelve-type source directory comes from 00_Config/paths.py.
CFG="$(cd "$MOD/../../00_Config" && pwd)"
SOUND12="$(python3 -c "import sys; sys.path.insert(0, '$CFG'); import paths; print(paths.SOUND_GSEA_DIR)")"
[ -d "$SOUND12" ] || { echo "merge_thirteen: paths.SOUND_GSEA_DIR is not a directory: $SOUND12" >&2; exit 1; }
OUT="$MOD/outputs/gsea_13types"

rm -rf "$OUT"; mkdir -p "$OUT"
cp "$SOUND12"/*_hallmark.csv "$OUT"/
n12=$(ls "$OUT" | wc -l)
cp "$MOD"/outputs/gsea/Neutrophils_*_hallmark.csv "$OUT"/
n13=$(ls "$OUT" | wc -l)
echo "$n12 Hallmark tables from the twelve-type run"
echo "$((n13 - n12)) added for neutrophils -> $n13 in $OUT"
