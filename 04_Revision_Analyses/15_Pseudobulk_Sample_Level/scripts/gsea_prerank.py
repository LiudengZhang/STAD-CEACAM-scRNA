"""
Prerank GSEA on the sample-level differential expression tables.

The metric, the gene set file and every GSEA setting are
12_R1.8_DEG_Recompute/scripts/recompute_deg.py's, taken from its `gsea()` at
lines 336-355, so that a NES here is comparable with the published per-cell
one:

    metric = logFC x -log10(P)          recompute_deg.py:343
    gene sets = 00_Reference/MSigDB_Hallmark_2020.gmt   (pinned, not Enrichr)
    permutation_num = 1000, min_size = 15, max_size = 500
    seed = 42                           the same SEED the published run used

Ties on gene symbol are resolved as there: the first row wins. Genes with no
p-value are dropped - DESeq2 sets them for Cook's outliers and independent
filtering - which is the same `dropna` the published run applies.

Every table is done in one process: 78 `conda run` startups on a loaded shared
node cost more than the enrichments do.

Shardable for the same reason run_deseq2.R is: <shard> of <nshards>, each
process taking every nshards-th table, and a table whose enrichment already
exists is skipped so an interrupted run resumes. The seed is per-table, so
sharding changes no result.

Run: python gsea_prerank.py <de_dir> <gsea_dir> [shard] [nshards]
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import HALLMARK_GMT                                       # noqa: E402

SEED = 42


def main():
    de, out = Path(sys.argv[1]), Path(sys.argv[2])
    shard = int(sys.argv[3]) if len(sys.argv) > 4 else 1
    nshards = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    out.mkdir(parents=True, exist_ok=True)
    import gseapy as gp
    print(f"gseapy {gp.__version__}, seed {SEED}, 1000 permutations, "
          f"{HALLMARK_GMT.name}", flush=True)

    tables = sorted(de.glob("*.csv"))
    for i, src in enumerate(tables):
        if i % nshards != (shard - 1) % nshards:
            continue
        base = src.stem                       # <cell>_<phase>_<method>
        dst = out / f"{base}_hallmark.csv"
        if dst.exists():
            print(f"{base}: have", flush=True)
            continue
        d = pd.read_csv(src).dropna(subset=["logfoldchanges", "pvals"]).copy()
        d["metric"] = d["logfoldchanges"] * -np.log10(d["pvals"].clip(lower=1e-300))
        rnk = (d.set_index("gene")["metric"].groupby(level=0).first()
               .sort_values(ascending=False))
        if len(rnk) < 15:
            print(f"{base}: only {len(rnk)} ranked genes - not enriched",
                  flush=True)
            continue
        res = gp.prerank(rnk=rnk, gene_sets=str(HALLMARK_GMT),
                         permutation_num=1000, min_size=15, max_size=500,
                         seed=SEED, no_plot=True, outdir=None, threads=2)
        r = res.res2d.copy()
        r["Term"] = r["Term"].astype(str).str.replace(r"^.*__", "", regex=True)
        r.to_csv(dst, index=False)
        print(f"{base:<34} {len(r):>3} Hallmark sets from {len(rnk):>6} "
              f"ranked genes", flush=True)
    print("gsea done")


if __name__ == "__main__":
    main()
