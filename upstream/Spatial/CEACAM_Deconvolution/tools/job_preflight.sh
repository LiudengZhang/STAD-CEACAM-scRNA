#!/bin/bash
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
# The filesystem roots below were specific to the machines this pipeline
# was run on. They are replaced by /path/to/machine and /path/to/home
# placeholders for deposition; repoint them at your own storage before
# running this script.
set -euo pipefail
ROOT=/path/to/home/graphst_repro
NTHREADS="$1"
echo "host=$(hostname) nproc=$(nproc) LSB_DJOB_NUMPROC=${LSB_DJOB_NUMPROC:-unset} requested=$NTHREADS"
[ "${LSB_DJOB_NUMPROC:-unset}" = "$NTHREADS" ] || { echo "FATAL: slot/thread mismatch"; exit 3; }
export OMP_NUM_THREADS=$NTHREADS OPENBLAS_NUM_THREADS=$NTHREADS MKL_NUM_THREADS=$NTHREADS
export NUMEXPR_NUM_THREADS=$NTHREADS VECLIB_MAXIMUM_THREADS=$NTHREADS
export PYTHONHASHSEED=42 GRAPHST_PIN=1 GRAPHST_NUM_THREADS=$NTHREADS GRAPHST_SEED=42 GRAPHST_DETERMINISTIC=1
export PYTHONPATH=$ROOT/tools/pysite PYTHONNOUSERSITE=1 MPLBACKEND=Agg
export MPLCONFIGDIR=$ROOT/run/.mplconfig; mkdir -p "$MPLCONFIGDIR"
"$ROOT/env/bin/python" -u "$ROOT/tools/preflight.py"
