#!/bin/bash
# Three stages, sequential. PYTHONHASHSEED is set here and nowhere else;
# gseapy is given an explicit seed in the script.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$(dirname "$HERE")/../.."
for stage in metrics settings seeds; do
  echo "########## $stage ##########"
  PYTHONHASHSEED=0 conda run --no-capture-output -n Liudeng_Python_310 \
      python "$HERE/metric_sensitivity.py" --stage "$stage" --threads 4
done
