#!/bin/bash
# The interpreter, R and library locations below were machine-specific
# paths under the author's home directory. They are replaced by a
# /path/to/home placeholder for deposition; point them at your own conda
# environment before running this script.
# drive_list.sh <spec> <jobfile> <streams> <threads_per_job>
# jobfile: one "<cell> <phase>" per line, in the order they should be started.
# Ordering matters for spec B: the cost of a glmer fit is roughly (cells x
# genes), and Epithelial and MoMac carry most of it, so they are started first
# and given the most threads.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
MOD="$(dirname "$HERE")"
ENV_LIB=/path/to/home/conda/envs/stad_ceacam/lib
SPEC="$1"; JOBS="$2"; STREAMS="$3"; THREADS="$4"
mkdir -p "$MOD/logs/$SPEC"
cat "$JOBS" | xargs -P "$STREAMS" -n 2 bash -c '
  c="$0"; p="$1"
  LD_LIBRARY_PATH='"$ENV_LIB"' OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  conda run --no-capture-output -n stad_ceacam python "'"$HERE"'/run_spec.py" \
    --cell "$c" --phase "$p" --spec "'"$SPEC"'" --threads "'"$THREADS"'" \
    > "'"$MOD"'/logs/'"$SPEC"'/${c}_${p}.log" 2>&1 \
  && echo "OK   $c $p" || echo "FAIL $c $p"
'
echo "=== $SPEC $(basename "$JOBS") finished ==="
