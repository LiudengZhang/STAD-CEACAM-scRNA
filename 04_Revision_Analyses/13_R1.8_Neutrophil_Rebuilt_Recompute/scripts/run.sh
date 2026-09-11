#!/bin/bash
# The interpreter, R and library locations below were machine-specific
# paths under the author's home directory. They are replaced by a
# /path/to/home placeholder for deposition; point them at your own conda
# environment before running this script.
# Neutrophils only, one phase per invocation, so the summary lands as
# recompute_summary_Neutrophils_<phase>.csv - the same file names the archived
# sound run produced for the other twelve types.
#
# LD_LIBRARY_PATH puts the environment's own libstdc++ ahead of the system one.
# Without it rpy2 cannot open libR.so, every MAST job is recorded unavailable,
# and the t-test half runs on regardless - a silent half-failure.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
ENV_LIB=/path/to/home/conda/envs/stad_ceacam/lib
for phase in pre post; do
  echo "=== $phase ==="
  LD_LIBRARY_PATH="$ENV_LIB" conda run --no-capture-output \
      -n stad_ceacam python "$HERE/recompute_deg_neutrophils.py" \
      --cell Neutrophils --phase "$phase"
done
