"""
The gate. For each of the 26 contrasts (13 cell types x pre/post) on BOTH
differential-expression sets, measure how much gene-set signal the contrast
carries at all, and how much of the TNFa/NF-kB number is permutation noise.

Everything except the permutation seed is held fixed at the published run's
settings: ranking metric logFC x -log10(P) (recompute_deg.py:343), the pinned
00_Reference/MSigDB_Hallmark_2020.gmt, permutation_num=1000, min_size=15,
max_size=500, PYTHONHASHSEED=0. Seed 42 is the published seed; 1, 7, 13, 101,
2026 and 777 are added to measure the noise. 26 x 2 x 7 = 364 prerank runs.

Writes outputs/signal_runs.csv (one row per run) and outputs/signal.csv (one
row per contrast per DE set).
"""
import argparse
import itertools
import sys
import time

import numpy as np
import pandas as pd
from gsea_common import (load_rank, run_one, CELLS, PHASES, NFKB_TERM, GMT,
                         SEEDS, PUB_SEED, OUT)


def measure(tab, rnk):
    r = tab[tab["Term"] == NFKB_TERM]
    nes_all = tab["NES"].astype(float)
    es_all = tab["ES"].astype(float)
    q_all = tab["FDR q-val"].astype(float)
    row = dict(n_genes_ranked=len(rnk), n_sets_tested=len(tab),
               max_abs_nes=float(nes_all.abs().max()),
               max_nes=float(nes_all.max()), min_nes=float(nes_all.min()),
               n_sets_q25=int((q_all < 0.25).sum()),
               n_sets_q05=int((q_all < 0.05).sum()),
               n_sets_nomp05=int((tab["NOM p-val"].astype(float) < 0.05).sum()),
               tnfa_overlap=len(GMT[NFKB_TERM] & set(rnk.index)))
    if len(r):
        r = r.iloc[0]
        es = float(r["ES"])
        row.update(nfkb_nes=float(r["NES"]), nfkb_es=es,
                   nfkb_es_sign="+" if es > 0 else "-",
                   nfkb_nom_p=float(r["NOM p-val"]),
                   nfkb_fdr_q=float(r["FDR q-val"]),
                   nfkb_rank=int(r["rank"]),
                   n_sets_same_es_sign=int((np.sign(es_all) == np.sign(es)).sum()),
                   nfkb_nes_degenerate=bool(abs(abs(float(r["NES"])) - 1.0) < 1e-9
                                            and float(r["NOM p-val"]) >= 1.0))
    else:
        row.update(nfkb_nes=np.nan, nfkb_es=np.nan, nfkb_es_sign="",
                   nfkb_nom_p=np.nan, nfkb_fdr_q=np.nan, nfkb_rank=np.nan,
                   n_sets_same_es_sign=np.nan, nfkb_nes_degenerate=False)
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    rows, t0 = [], time.time()
    for degset in ("live", "sound13"):
        for cell, phase in itertools.product(CELLS, PHASES):
            rnk = load_rank(cell, phase, degset)
            for seed in SEEDS:
                tab = run_one(rnk, seed, threads=args.threads)
                rows.append(dict(degset=degset, cell=cell, phase=phase,
                                 seed=seed, **measure(tab, rnk)))
                print(f"[{time.time()-t0:6.0f}s] {degset:<8}{cell:<20}{phase:<5}"
                      f"s{seed:<5} NES={rows[-1]['nfkb_nes']:.4f} "
                      f"rank={rows[-1]['nfkb_rank']}", flush=True)
    runs = pd.DataFrame(rows)
    runs.to_csv(OUT / "signal_runs.csv", index=False)
    print(f"done in {(time.time()-t0)/60:.1f} min", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
