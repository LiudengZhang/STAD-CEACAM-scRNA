#!/bin/bash
# The same three stages on the `live` DEG tables (12_R1.8_DEG_Recompute/outputs/deg),
# which is the set the main text's counted claims are counted over. Waits for the
# sound13 sweep to finish so the two do not contend for CPU.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
while pgrep -f "metric_sensitivity.py --stage" > /dev/null; do sleep 20; done
for stage in metrics settings seeds; do
  echo "########## live $stage ##########"
  PYTHONHASHSEED=0 conda run --no-capture-output -n stad_ceacam \
      python "$HERE/metric_sensitivity.py" --stage "$stage" --degset live --threads 4
done
