#!/bin/bash
# limma-voom, then DESeq2, then prerank GSEA on both, for every contrast
# build_pseudobulk.py wrote.
#
# Two R environments, because no single one here has both halves: limma 3.62.2
# and edgeR 4.4.0 are in r_demo, DESeq2 1.50.2 is in deseq2_env, and neither has
# the other. Both were checked before this script was written; a missing package
# fails loudly.
#
# One session per stage rather than one per contrast. The first version called
# `conda run` 78 times and spent nine minutes on the first contrast.
#
# The two slow stages are sharded four ways. The node was shared with a job
# holding 26 of its 48 cores and DESeq2 was running at about five minutes a
# contrast; four is a share of what was left, not a claim on the machine. Each
# stage skips a contrast whose output already exists, so a killed run resumes.
# Sharding changes no result: the contrasts are independent fits and the GSEA
# seed is per table.
#
# limma runs first because it writes the filterByExpr gene list DESeq2 reads.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
PB="$HERE/../outputs/pseudobulk"
DE="$HERE/../outputs/de"
GS="$HERE/../outputs/gsea"
LOG="$HERE/../logs"
mkdir -p "$DE" "$GS" "$LOG"
N=4

conda run --no-capture-output -n r_demo Rscript "$HERE/run_limma_voom.R" "$PB" "$DE"

for i in $(seq 1 $N); do
  conda run --no-capture-output -n deseq2_env Rscript "$HERE/run_deseq2.R" \
      "$PB" "$DE" "$i" "$N" > "$LOG/deseq2_shard$i.log" 2>&1 &
done
wait

for i in $(seq 1 $N); do
  conda run --no-capture-output -n Liudeng_Python_310 python \
      "$HERE/gsea_prerank.py" "$DE" "$GS" "$i" "$N" \
      > "$LOG/gsea_shard$i.log" 2>&1 &
done
wait

echo "=== all contrasts done"
