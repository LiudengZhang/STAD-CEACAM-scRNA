#!/usr/bin/env bash
# The interpreter, R and library locations below were machine-specific
# paths under the author's home directory. They are replaced by a
# /path/to/home placeholder for deposition; point them at your own conda
# environment before running this script.
# Launch all 26 recomputations, one process per (cell type, timepoint).
# Separate processes because MAST runs through rpy2 and R is not thread-safe.
# Each job needs about 13 GB of memory; all 26 start at once.
#   bash work/launch_all.sh               run everything
#   bash work/launch_all.sh --self-test   test the failure reporting only
set -u
cd "$(dirname "$0")/.."

# Wait on every launched job by PID (a bare `wait` returns 0 regardless).
# Takes "pid:label" words, prints the failures, returns 1 if there were any.
wait_all() {
    local failed= n=0
    for job in "$@"; do
        local pid=${job%%:*} what=${job#*:} rc=0
        wait "$pid" || rc=$?
        if [ "$rc" -ne 0 ]; then
            failed="$failed\n  $what  rc=$rc"
            n=$((n + 1))
        fi
    done
    if [ "$n" -gt 0 ]; then
        printf 'NOT ALL DONE: %d of %d did not succeed:%b\n' \
               "$n" "$#" "$failed" >&2
        echo "their logs are under logs/." >&2
        return 1
    fi
    return 0
}

# Self-test: two jobs, one of which fails, must be reported as a failure.
if [ "${1:-}" = "--self-test" ]; then
    st=0
    ( exit 0 ) & a=$!
    ( exit 3 ) & b=$!
    if wait_all "$a:good_job" "$b:bad_job" 2>/tmp/wait_all.$$; then
        echo "FAIL: one job exited 3 and wait_all reported success"; st=1
    else
        grep -q 'bad_job  rc=3' /tmp/wait_all.$$ \
            || { echo "FAIL: the failure was not named with its status"; st=1; }
        grep -q '1 of 2' /tmp/wait_all.$$ \
            || { echo "FAIL: the count was not reported"; st=1; }
    fi
    ( exit 0 ) & c=$!
    wait_all "$c:only_good" || { echo "FAIL: a clean run was called failed"; st=1; }
    rm -f /tmp/wait_all.$$
    [ "$st" -eq 0 ] && echo "wait_all behaves as required" \
                    || echo "WAIT_ALL DID NOT BEHAVE AS REQUIRED" >&2
    exit "$st"
fi

PY=/path/to/home/conda/envs/stad_ceacam/bin/python
export PYTHONNOUSERSITE=1
export R_HOME=/path/to/home/conda/envs/r_bayesprism/lib/R
export LD_LIBRARY_PATH=/path/to/home/conda/envs/r_bayesprism/lib:/path/to/home/conda/envs/r_bayesprism/lib/R/lib
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2

CELLS="B_cells DC_cells Endothelial_cells Epithelial Fibroblast Mast_cells MoMac
       Neutrophils NK_cells Pericyte Plasma_cells TCD4_cells TCD8_cells"
JOBS=()
for c in $CELLS; do
  for p in pre post; do
    nohup $PY work/recompute_deg.py --cell "$c" --phase "$p" \
      > "logs/${c}_${p}.log" 2>&1 &
    echo "launched $c $p  pid $!"
    JOBS+=("$!:$c $p")
  done
done
echo "--- $(jobs -rp | wc -l) processes running ---"
wait_all "${JOBS[@]}" || exit 1
echo "ALL DONE $(date)"
