"""
The paragraph-70 / Figure 5N-5I numbers, which no sibling module has touched and
which verify_numbers.py does not check.

    Epithelial cells showed enrichment for inflammatory response (NES = 1.48,
    P = 0.02) and hypoxia (NES = 1.56, P = 0.006) ... Fibroblasts exhibited
    enrichment for inflammatory response (NES = 1.30, P = 0.02) ...
    Monocytes/macrophages ... inflammatory response (NES = 1.67), EMT
    (NES = 1.66), hypoxia (NES = 1.53), and angiogenesis (NES = 1.50; all
    P <= 0.05)

These come from a DIFFERENT pipeline from every other GSEA number in the paper:
the differential expression is computed inside the panel scripts
(05_G/create_panel_g_gsea_2types.py, 05_GSEA_Summary/run_gsea_momac.py,
run_gsea_5celltypes.py), with no 10%-detection gene filter, and the prerank uses
**min_size = 5**, not the min_size = 15 of the NF-kB analysis.

This script copies those scripts' own DE and prerank calls exactly - same obs
columns, same 5,000-cell subsample at random_state 42, same rank metric, same
GMT, same settings - and repeats the prerank at seven seeds. Nothing is written
outside this module; the panels' own directories are read, never touched.

`rank_genes_groups` is called without `use_raw`, exactly as the panel scripts
call it; all three objects carry `.raw`, so scanpy's default reads `.raw` and
not the doubly-normalised `.X`. That is checked and recorded here.
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import (MOMAC_H5AD, EPITHELIAL_H5AD, FIBROBLAST_H5AD,
                   DC_CELLS_H5AD, HALLMARK_GMT)   # noqa: E402
from gsea_common import OUT, SEEDS                # noqa: E402

COMPARISON = dict(column="stomach_post_grouping", group1="Responsed",
                  group2="No-response")
TARGETS = {
    ("MoMac", "Inflammatory Response"): 1.67,
    ("MoMac", "Epithelial Mesenchymal Transition"): 1.66,
    ("MoMac", "Hypoxia"): 1.53,
    ("MoMac", "Angiogenesis"): 1.50,
    ("Epithelial", "Inflammatory Response"): 1.48,
    ("Epithelial", "Hypoxia"): 1.56,
    ("Fibroblast", "Inflammatory Response"): 1.30,
}
SOURCES = {"MoMac": MOMAC_H5AD, "Epithelial": EPITHELIAL_H5AD,
           "Fibroblast": FIBROBLAST_H5AD, "DC": DC_CELLS_H5AD}

rows, allterms = [], []
for cell, path in SOURCES.items():
    print(f"=== {cell} ===", flush=True)
    ad = sc.read_h5ad(path)
    ad.obs_names_make_unique()
    used_raw = ad.raw is not None
    ad = ad[ad.obs["Sample site"] == "Stomach"].copy()
    ad = ad[ad.obs["Treatment phase"] == "Post"].copy()
    col, g1, g2 = COMPARISON["column"], COMPARISON["group1"], COMPARISON["group2"]
    sub = ad[ad.obs[col].isin([g1, g2])].copy()
    n1, n2 = int((sub.obs[col] == g1).sum()), int((sub.obs[col] == g2).sum())
    if sub.n_obs > 5000:
        sc.pp.subsample(sub, n_obs=5000, random_state=42)
    sc.tl.rank_genes_groups(sub, groupby=col, groups=[g2], reference=g1,
                            method="t-test_overestim_var", pts=True)
    d = sc.get.rank_genes_groups_df(sub, group=g2)
    d["rank_metric"] = d["logfoldchanges"] * -np.log10(d["pvals"].clip(lower=1e-300))
    rnk = d.set_index("names")["rank_metric"].dropna()
    rnk = rnk[~rnk.index.duplicated(keep="first")].sort_values(ascending=False)
    print(f"    R {n1} NR {n2}, ranked list {len(rnk)}, .raw used = {used_raw}",
          flush=True)
    for seed in SEEDS:
        res = gp.prerank(rnk=rnk, gene_sets=str(HALLMARK_GMT), min_size=5,
                         max_size=500, permutation_num=1000, seed=seed,
                         no_plot=True, outdir=None, threads=4)
        t = res.res2d.copy()
        t["Term"] = t["Term"].astype(str).str.replace(r"^.*__", "", regex=True)
        for c in ("ES", "NES", "NOM p-val", "FDR q-val"):
            t[c] = pd.to_numeric(t[c], errors="coerce")
        t = t.sort_values("NES", ascending=False).reset_index(drop=True)
        t["rank"] = np.arange(1, len(t) + 1)
        t["cell"], t["seed"] = cell, seed
        allterms.append(t[["cell", "seed", "Term", "ES", "NES", "NOM p-val",
                           "FDR q-val", "rank"]])
        for (c2, pw), quoted in TARGETS.items():
            if c2 != cell:
                continue
            r = t[t["Term"] == pw]
            if r.empty:
                continue
            r = r.iloc[0]
            rows.append(dict(cell=cell, pathway=pw, quoted_nes=quoted, seed=seed,
                             measured_nes=round(float(r["NES"]), 4),
                             nom_p=round(float(r["NOM p-val"]), 4),
                             fdr_q=round(float(r["FDR q-val"]), 4),
                             rank=int(r["rank"]), n_sets=len(t),
                             n_genes_ranked=len(rnk), raw_used=used_raw,
                             n_cells_R=n1, n_cells_NR=n2))
        print(f"    seed {seed} done", flush=True)
    del ad, sub

d = pd.DataFrame(rows)
d.to_csv(OUT / "para70_reproduction.csv", index=False)
pd.concat(allterms, ignore_index=True).to_csv(OUT / "para70_all_terms.csv", index=False)

s = (d.groupby(["cell", "pathway", "quoted_nes"])
       .agg(nes_seed42=("measured_nes", lambda x: x.iloc[0]),
            nes_min=("measured_nes", "min"), nes_max=("measured_nes", "max"),
            p_seed42=("nom_p", lambda x: x.iloc[0]),
            q_min=("fdr_q", "min"), q_max=("fdr_q", "max"),
            rank_min=("rank", "min"), rank_max=("rank", "max"),
            n_sets=("n_sets", "max")).reset_index())
s["abs_diff_from_printed"] = (s["nes_seed42"] - s["quoted_nes"]).abs().round(4)
s.to_csv(OUT / "para70_summary.csv", index=False)
pd.set_option("display.width", 220)
print(s.to_string(index=False))
sys.exit(0)
