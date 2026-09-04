#!/bin/bash
# The filesystem roots below were specific to the machines this pipeline
# was run on. They are replaced by /path/to/machine and /path/to/home
# placeholders for deposition; repoint them at your own storage before
# running this script.
# =============================================================================
# GraphST spatial deconvolution - reproduction driver (CIR-26-0753-ET, WP3)
#
# This file is the record of EVERY pinned parameter.  The staged analysis
# scripts under code/ are a byte-identical provenance copy and are never
# edited; everything below is imposed from outside them, through
#   - the rewritten config copies      (paths + explicit device only)
#   - tools/pysite/sitecustomize.py    (threads, seeds, deterministic algs)
#   - the environment variables set in tools/job_graphst.sh
#
# Lane: Seadragon LSF, CPU.  The original run was CPU (its log records
# "No GPU detected - using CPU"), so no GPU is requested.
#
#   ./run.sh smoke    submit the two byte-equality smoke jobs (sample_08)
#   ./run.sh compare  byte-compare the two smoke outputs
#   ./run.sh fanout   submit the remaining 9 samples as an idempotent array
#   ./run.sh params   print the pinned parameters
# =============================================================================
set -euo pipefail
ROOT=/path/to/home/graphst_repro

# ---------------------------------------------------------------- pinned ----
NTHREADS=24          # LSF -n and OMP/OPENBLAS/MKL/NUMEXPR/VECLIB and
                     # torch.set_num_threads / set_num_interop_threads all == this.
                     # 24 = torch.get_num_threads() default on the original
                     # 48-logical-core host, i.e. the closest available match to
                     # the original run's (never recorded) thread count.
SEED=42              # config random_seed, PYTHONHASHSEED, random/numpy/torch.
DEVICE=cpu           # set explicitly in the config copies; detect_device("auto")
                     # is deliberately NOT relied on.
DETERMINISTIC=1      # torch.use_deterministic_algorithms(True)

SMOKE_SAMPLE=sample_08     # fastest sample in the original run (123.0 min)

# sized from the smoke measurement -- see VERDICT.md
HOSTGRP=comlow       # 28-core/188G Intel_Platinum pool; keeps smoke_a and smoke_b
                     # on identical hardware so a byte difference cannot be a
                     # CPU-dispatch artefact, and keeps the fan-out comparable.
QUEUE=${QUEUE:-short}  # smoke only.  prio 76, 180 min max wall, backlog of 5 when
                     # probed; the smoke measured 3 s of queue wait.
                     # NEVER 'medium' (~2.9M queued).
WALL=${WALL:-175}    # minutes, -W, for the smoke
MEM=150              # GB, -M / rusage[mem]

# Fan-out routing, decided by MEASUREMENT rather than by preference.  On
# Seadragon (24 threads, comlow) the smoke reproduced the original host's
# training rate to within 5% -- phase 1 25.3 vs 26.6 it/s, phase 2 2.20 vs
# 2.15 s/it -- so the original per-sample wall times transfer directly:
#   01 347  02 331  03 198  04 242  05 459  06 224  07 160  08 123  09 193  10 191  (min)
# Only sample_08 (123) and sample_07 (160) are under the 180 min 'short'
# ceiling, and 160 leaves 9% of headroom against a 5% rate uncertainty.  A
# wall-clock kill costs the whole 3 h, so the remaining nine all go to 'long'.
QUEUE_FANOUT=long    # prio 25.  Backlog is large but it drains (7797 pending
                     # against 7737 running when probed).
WALL_FANOUT=1500     # minutes.  NOT a size estimate: 'long' has a wall-limit
                     # FLOOR, not just a ceiling -- its esub rejects anything
                     # with RUNLIMIT <= 86400 s (1440 min), verbatim
                     #   "RUNLIMIT must be greater than 86400 seconds."
                     # Measured bands: short <=180 min, medium <=1440 min,
                     # long >1440..14400 min, vlong <=30240 min.  A 2-8 h job
                     # therefore has no home except 'medium' (~2.9M queued,
                     # excluded) so it goes to 'long' at the smallest wall the
                     # queue will accept.  Worst case here is ~1000 min: the
                     # slowest original sample (459 min) on a node like the one
                     # the first smoke_b landed on, which ran 2.2x slower than
                     # the original host for no visible reason (its load was
                     # LOWER than the fast node's).

# unchanged analysis parameters, restated here so a diff is visible:
#   epochs 1200 | learning_rate 0.001 | deconvolution true | n_neighbors 6
#   retain_percent 0.15 | clustering leiden | resolution 0.5
#   n_clusters_range [6..12] | reference stomach_14types_reference.h5ad
#   cell_type_column cell_type_14
# -----------------------------------------------------------------------------

params() {
    sed -n '/^# ------.*pinned ----$/,/^# ----------------------------$/p' "$0" | sed 's/^/  /'
    echo "  env prefix : $ROOT/env"
    echo "  code md5s  : $ROOT/manifests/code_md5.txt"
}

submit_one() {   # submit_one <run_name> <sample> <jobname>
    local run="$1" sample="$2" name="$3"
    mkdir -p "$ROOT/logs/lsf" "$ROOT/run/$run/results"
    bsub -q "$QUEUE" -m "$HOSTGRP" -n "$NTHREADS" -M "$MEM" \
         -R "rusage[mem=$MEM] span[hosts=1]" -W "$WALL" -J "$name" \
         -o "$ROOT/logs/lsf/${name}.out" -e "$ROOT/logs/lsf/${name}.err" \
         "$ROOT/tools/job_graphst.sh" "$run" "$sample" "$NTHREADS"
}

case "${1:-params}" in
  params) params ;;
  smoke)
    submit_one smoke_a "$SMOKE_SAMPLE" gstA
    submit_one smoke_b "$SMOKE_SAMPLE" gstB
    ;;
  smoke_a) submit_one smoke_a "$SMOKE_SAMPLE" gstA ;;
  smoke_b) submit_one smoke_b "$SMOKE_SAMPLE" gstB ;;
  compare)
    A="$ROOT/run/smoke_a/results/${SMOKE_SAMPLE}_graphst_output.h5ad"
    B="$ROOT/run/smoke_b/results/${SMOKE_SAMPLE}_graphst_output.h5ad"
    ls -l "$A" "$B"
    if cmp -s "$A" "$B"; then echo "BYTE_IDENTICAL yes"; else echo "BYTE_IDENTICAL no"; cmp "$A" "$B" || true; fi
    md5sum "$A" "$B"
    ;;
  fanout)
    # 9 remaining samples, idempotent (output exists => exit 0), so a resubmit
    # is a resume.  %3 throttles concurrency; no automatic requeue, so an
    # unknown exit code fails closed.
    mkdir -p "$ROOT/logs/lsf" "$ROOT/run/full/results"
    bsub -q "$QUEUE_FANOUT" -m "$HOSTGRP" -n "$NTHREADS" -M "$MEM" \
         -R "rusage[mem=$MEM] span[hosts=1]" -W "$WALL_FANOUT" -J "gstF[1-10]%3" \
         -o "$ROOT/logs/lsf/gstF.%I.out" -e "$ROOT/logs/lsf/gstF.%I.err" \
         "$ROOT/tools/job_graphst.sh" full "@INDEX" "$NTHREADS"
    ;;
  *) echo "usage: [QUEUE=q] [WALL=min] $0 {params|smoke|smoke_a|smoke_b|compare|fanout}"; exit 2 ;;
esac
