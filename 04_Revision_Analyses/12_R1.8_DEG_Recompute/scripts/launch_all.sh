#!/usr/bin/env bash
# The interpreter, R and library locations below were machine-specific
# paths under the author's home directory. They are replaced by a
# /path/to/home placeholder for deposition; point them at your own conda
# environment before running this script.
# Launch all 26 recomputations, one process per (cell type, timepoint).
#
# One process per job rather than a thread pool: MAST goes through rpy2 into a
# single R interpreter per process, and R is not thread-safe. 12.6 GB each,
# 26 at once is 330 GB against the 685 GB free, so they all start at once.
set -u
cd "$(dirname "$0")/.."
PY=/path/to/home/miniforge3/envs/stad_final2/bin/python
export PYTHONNOUSERSITE=1
export R_HOME=/path/to/home/miniforge3/envs/r_demo/lib/R
export LD_LIBRARY_PATH=/path/to/home/miniforge3/envs/r_demo/lib:/path/to/home/miniforge3/envs/r_demo/lib/R/lib
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2

CELLS="B_cells DC_cells Endothelial_cells Epithelial Fibroblast Mast_cells MoMac
       Neutrophils NK_cells Pericyte Plasma_cells TCD4_cells TCD8_cells"
for c in $CELLS; do
  for p in pre post; do
    nohup $PY scripts/recompute_deg.py --cell "$c" --phase "$p" \
      > "logs/${c}_${p}.log" 2>&1 &
    echo "launched $c $p  pid $!"
  done
done
echo "--- $(jobs -rp | wc -l) processes running ---"
wait
echo "ALL DONE $(date)"
