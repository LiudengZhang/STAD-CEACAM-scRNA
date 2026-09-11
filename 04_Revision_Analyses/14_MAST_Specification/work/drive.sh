#!/bin/bash
# The interpreter, R and library locations below were machine-specific
# paths under the author's home directory. They are replaced by a
# /path/to/home placeholder for deposition; point them at your own conda
# environment before running this script.
# One process per (cell, phase, spec). Sequential inside a stream, several
# streams in parallel; the stream count and the per-job R thread count are
# arguments because this machine is shared with three other agents.
#
# LD_LIBRARY_PATH puts the environment's own libstdc++ ahead of the system
# one. Without it rpy2 cannot open libR.so and MAST is silently unavailable.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
MOD="$(dirname "$HERE")"
ENV_LIB=/path/to/home/conda/envs/stad_ceacam/lib
SPEC="$1"; STREAMS="${2:-6}"; THREADS="${3:-1}"
CELLS="B_cells DC_cells Endothelial_cells Epithelial Fibroblast Mast_cells MoMac Neutrophils NK_cells Pericyte Plasma_cells TCD4_cells TCD8_cells"
mkdir -p "$MOD/logs/$SPEC"
for c in $CELLS; do for p in pre post; do echo "$c $p"; done; done | \
xargs -P "$STREAMS" -n 2 bash -c '
  c="$0"; p="$1"
  LD_LIBRARY_PATH='"$ENV_LIB"' OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  conda run --no-capture-output -n stad_ceacam python "'"$HERE"'/run_spec.py" \
    --cell "$c" --phase "$p" --spec "'"$SPEC"'" --threads "'"$THREADS"'" \
    > "'"$MOD"'/logs/'"$SPEC"'/${c}_${p}.log" 2>&1 \
  && echo "OK   $c $p" || echo "FAIL $c $p"
'
echo "=== $SPEC driver finished ==="
