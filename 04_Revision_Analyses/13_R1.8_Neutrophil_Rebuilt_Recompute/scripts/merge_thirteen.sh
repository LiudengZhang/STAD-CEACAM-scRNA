#!/bin/bash
# The thirteen-type sound table: the twelve types of the 2026-08-31 sound run,
# plus neutrophils from the rebuilt input. Nothing is merged in place - the
# archive is copied from, never written to, and 12_R1.8_DEG_Recompute/outputs/
# is not touched.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
MOD="$(dirname "$HERE")"
ROOT="$(cd "$MOD/../.." && pwd)"
SOUND12="$ROOT/07_Archive/2026-08-31_deg_recompute_on_sound_per_cell_type_inputs/02_New_Analyses/12_R1.8_DEG_Recompute/outputs/gsea"
OUT="$MOD/outputs/gsea_13types"

rm -rf "$OUT"; mkdir -p "$OUT"
cp "$SOUND12"/*_hallmark.csv "$OUT"/
n12=$(ls "$OUT" | wc -l)
cp "$MOD"/outputs/gsea/Neutrophils_*_hallmark.csv "$OUT"/
n13=$(ls "$OUT" | wc -l)
echo "$n12 Hallmark tables from the twelve-type sound run"
echo "$((n13 - n12)) added for neutrophils -> $n13 in $OUT"
