#!/bin/bash
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
# The filesystem roots below were specific to the machines this pipeline
# was run on. They are replaced by /path/to/machine and /path/to/home
# placeholders for deposition; repoint them at your own storage before
# running this script.
# GraphST reproduction - one sample, one LSF job.
# usage: job_graphst.sh <run_name> <sample_key> <nthreads>
set -euo pipefail

ROOT=/path/to/home/graphst_repro
RUN="$1"; SAMPLE="$2"; NTHREADS="$3"

# In an LSF array the sample key comes from the index.  Resolving it here
# (rather than in the bsub line) keeps the submission free of shell quoting
# that bsub would pass through literally.
if [ "$SAMPLE" = "@INDEX" ]; then
    [ -n "${LSB_JOBINDEX:-}" ] || { echo "FATAL: @INDEX outside a job array"; exit 5; }
    SAMPLE=$(printf "sample_%02d" "$LSB_JOBINDEX")
fi
echo "RUN=$RUN SAMPLE=$SAMPLE"

OUTDIR="$ROOT/run/$RUN/results"
OUT="$OUTDIR/${SAMPLE}_graphst_output.h5ad"
mkdir -p "$OUTDIR"

# --- idempotent: an existing output means this index is done -----------------
if [ -s "$OUT" ]; then
    echo "IDEMPOTENT_SKIP $OUT ($(stat -c %s "$OUT") bytes)"
    exit 0
fi

# --- never trust nproc: assert LSF actually gave us what we pinned to --------
echo "host=$(hostname) nproc=$(nproc) LSB_DJOB_NUMPROC=${LSB_DJOB_NUMPROC:-unset} requested=$NTHREADS"
if [ "${LSB_DJOB_NUMPROC:-unset}" != "$NTHREADS" ]; then
    echo "FATAL: LSF allocated ${LSB_DJOB_NUMPROC:-unset} slots but threads are pinned to $NTHREADS"
    exit 3
fi

# --- pinned determinism (also echoed by tools/pysite/sitecustomize.py) -------
export OMP_NUM_THREADS="$NTHREADS"
export OPENBLAS_NUM_THREADS="$NTHREADS"
export MKL_NUM_THREADS="$NTHREADS"
export NUMEXPR_NUM_THREADS="$NTHREADS"
export VECLIB_MAXIMUM_THREADS="$NTHREADS"
export PYTHONHASHSEED=42
export GRAPHST_PIN=1
export GRAPHST_NUM_THREADS="$NTHREADS"
export GRAPHST_SEED=42
export GRAPHST_DETERMINISTIC=1
export PYTHONPATH="$ROOT/tools/pysite"
export PYTHONNOUSERSITE=1
export MPLBACKEND=Agg
export MPLCONFIGDIR="$ROOT/run/$RUN/.mplconfig"
mkdir -p "$MPLCONFIGDIR"

cd "$ROOT/code/02_GraphST_Analysis/scripts"

T0=$(date +%s)
"$ROOT/env/bin/python" -u 01_run_graphst_per_sample.py \
    "$ROOT/run/$RUN/config/graphst_config.yaml" "$SAMPLE" &
PY=$!

# peak RSS from VmHWM of the python process itself.  ru_maxrss of a child can be
# the PARENT's high-water mark (see the compute-routing TRAPS table), so read
# /proc/<pid>/status directly.
PEAK=0
while kill -0 "$PY" 2>/dev/null; do
    V=$(awk '/VmHWM/{print $2}' /proc/$PY/status 2>/dev/null || echo 0)
    [ -n "$V" ] && [ "$V" -gt "$PEAK" ] 2>/dev/null && PEAK=$V
    sleep 15
done
wait "$PY"; RC=$?
T1=$(date +%s)

echo "PEAK_VmHWM_KB=$PEAK  PEAK_GiB=$(awk -v k=$PEAK 'BEGIN{printf "%.2f", k/1048576}')"
echo "WALL_SECONDS=$((T1-T0))  WALL_MINUTES=$(awk -v s=$((T1-T0)) 'BEGIN{printf "%.1f", s/60}')"
echo "PYTHON_EXIT=$RC"
[ "$RC" -eq 0 ] || exit "$RC"
[ -s "$OUT" ] || { echo "FATAL: python exited 0 but $OUT is missing/empty"; exit 4; }
echo "OUTPUT_BYTES=$(stat -c %s "$OUT")"
md5sum "$OUT"
echo "JOB_OK $RUN $SAMPLE"
