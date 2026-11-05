#!/bin/bash
# The interpreter, R and library locations below were machine-specific
# paths under the author's home directory. They are replaced by a
# /path/to/home placeholder for deposition; point them at your own conda
# environment before running this script.
# Neutrophil DE and GSEA, one phase per call
# (writes recompute_summary_Neutrophils_<phase>.csv). Run: bash run.sh
#
# ENV_LIB must point at your conda environment's lib directory: it puts that
# libstdc++ ahead of the system one so rpy2 can load R for MAST. Without it
# MAST is silently skipped, so the path is checked first.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
ENV_LIB=/path/to/home/conda/envs/stad_ceacam/lib
if [ ! -d "$ENV_LIB" ]; then
    echo "ENV_LIB does not exist: $ENV_LIB" >&2
    echo "Point it at the lib directory of the conda environment you are" >&2
    echo "running in. Refusing to continue: an absent LD_LIBRARY_PATH is" >&2
    echo "ignored silently and MAST would be recorded unavailable for every" >&2
    echo "job while the t-test half ran on and the run exited 0." >&2
    exit 2
fi
if [ ! -e "$ENV_LIB/libstdc++.so.6" ]; then
    echo "$ENV_LIB exists but holds no libstdc++.so.6" >&2
    echo "That is the one library this path is here to put ahead of the" >&2
    echo "system's, so pointing LD_LIBRARY_PATH at it would achieve nothing." >&2
    exit 2
fi
for phase in pre post; do
  echo "=== $phase ==="
  LD_LIBRARY_PATH="$ENV_LIB" conda run --no-capture-output \
      -n stad_ceacam python "$HERE/recompute_deg_neutrophils.py" \
      --cell Neutrophils --phase "$phase"
done
