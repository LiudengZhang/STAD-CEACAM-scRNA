#!/bin/bash
# Rebuild every restyled S7-S11 panel. Version A is never written to.
set -u
cd "$(dirname "$0")"
R="conda run -n Liudeng_Python_310 python"
echo "=== standalone panel scripts (5) ==="
for f in S7_Cohort_Statistics/S7_A/create_S7_A_cohort_design.py \
         S7_Cohort_Statistics/S7_B/create_S7_B_forest.py \
         S8_CEACAM_Metaprogram/S8_E/create_S8_E_convergence.py \
         S8_CEACAM_Metaprogram/S8_F/create_S8_F_leave_one_out.py \
         S8_CEACAM_Metaprogram/S8_G/create_S8_G_rna_protein.py ; do
  [ -f "$f" ] || { echo "  MISSING $f"; continue; }
  $R "$f" 2>&1 | grep -viE 'futurewarning|pynvml|^ *import pynvml' | grep -E 'Saved|WARNING|Error|Traceback' | sed 's/^/  /'
done
echo
echo "=== drivers (analysis-drawn panels) ==="
cd _drivers
for f in restyle_*.py ; do
  echo "  -- $f"
  $R "$f" 2>&1 | grep -viE 'futurewarning|pynvml|^ *import pynvml' | grep -E 'mm$|WARNING|Error|Traceback|CONTENT' | sed 's/^/    /'
done
