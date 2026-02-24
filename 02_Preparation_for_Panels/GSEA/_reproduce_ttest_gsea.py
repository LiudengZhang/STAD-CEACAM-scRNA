#!/usr/bin/env python3
"""
Reproduce t-test GSEA radar data exactly.
Method: inline t-test DEGs (NR vs R), logFC × -log10(p), 100 permutations, min_size=5.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import gseapy as gp
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_Config"))
from paths import *

OUTPUT_DIR = Path(__file__).parent
PATHWAY_NAME = "TNF-alpha Signaling via NF-kB"
MAX_PATHWAYS = 50

CELL_TYPE_CONFIG = {
    "B_cells":           {"path": B_CELLS_H5AD,        "label": "B cells"},
    "DC_cells":          {"path": DC_CELLS_H5AD,        "label": "DC"},
    "Endothelial_cells": {"path": ENDOTHELIAL_H5AD,     "label": "Endothelial"},
    "Epithelial":        {"path": EPITHELIAL_H5AD,       "label": "Epithelial"},
    "Fibroblast":        {"path": FIBROBLAST_H5AD,       "label": "Fibroblast"},
    "Mast_cells":        {"path": MAST_CELLS_H5AD,      "label": "Mast"},
    "MoMac":             {"path": MOMAC_H5AD,            "label": "MoMac"},
    "Neutrophils":       {"path": NEUTROPHILS_H5AD,      "label": "Neutrophils"},
    "NK_cells":          {"path": NK_CELLS_H5AD,         "label": "NK"},
    "Pericyte":          {"path": PERICYTE_H5AD,         "label": "Pericyte"},
    "Plasma_cells":      {"path": PLASMA_CELLS_H5AD,    "label": "Plasma"},
    "TCD4_cells":        {"path": TCD4_H5AD,             "label": "CD4+ T"},
    "TCD8_cells":        {"path": TCD8_H5AD,             "label": "CD8+ T"},
}

# NR vs R
COMPARISONS = {
    "pre":  {"column": "stomach_pre_grouping",  "group1": "Responsed", "group2": "No-response"},
    "post": {"column": "stomach_post_grouping", "group1": "Responsed", "group2": "No-response"},
}


def run_deg_ttest(adata, comp_info, max_cells=5000):
    col = comp_info["column"]
    g1, g2 = comp_info["group1"], comp_info["group2"]
    if col not in adata.obs.columns:
        return None
    valid = adata.obs[col].isin([g1, g2])
    sub = adata[valid].copy()
    n1 = (sub.obs[col] == g1).sum()
    n2 = (sub.obs[col] == g2).sum()
    print(f"    R={n1}, NR={n2}")
    if n1 < 10 or n2 < 10:
        return None
    if sub.n_obs > max_cells:
        sc.pp.subsample(sub, n_obs=max_cells, random_state=42)
    try:
        sc.tl.rank_genes_groups(sub, groupby=col, groups=[g2], reference=g1,
                                method='t-test_overestim_var', pts=True)
        return sc.get.rank_genes_groups_df(sub, group=g2)
    except Exception as e:
        print(f"    DEG error: {e}")
        return None


def run_gsea(deg_df):
    if deg_df is None or len(deg_df) == 0:
        return None
    deg_df = deg_df.copy()
    deg_df['rank_metric'] = deg_df['logfoldchanges'] * (-np.log10(deg_df['pvals'].clip(lower=1e-300)))
    rnk = deg_df.set_index('names')['rank_metric'].dropna()
    rnk = rnk[~rnk.index.duplicated(keep='first')].sort_values(ascending=False)
    if len(rnk) < 15:
        return None
    try:
        res = gp.prerank(rnk=rnk, gene_sets='MSigDB_Hallmark_2020',
                         outdir=None, min_size=5, max_size=500,
                         permutation_num=100, seed=42, verbose=False)
        return res.res2d
    except Exception as e:
        print(f"    GSEA error: {e}")
        return None


def main():
    print("=" * 60)
    print("Reproducing t-test GSEA (NR vs R, 100 perms)")
    print("=" * 60)

    results = []
    for ct, cfg in CELL_TYPE_CONFIG.items():
        print(f"\n=== {ct} ===")
        row = {"cell_type": ct, "label": cfg["label"]}
        h5ad_path = Path(cfg["path"])
        if not h5ad_path.exists():
            print(f"  NOT FOUND: {h5ad_path}")
            for c in COMPARISONS:
                row[f"{c}_rank"] = None; row[f"{c}_nes"] = None; row[f"{c}_inverted"] = None
            results.append(row)
            continue

        adata = sc.read_h5ad(h5ad_path)
        adata.obs_names_make_unique()
        print(f"  {adata.n_obs} cells")

        for comp_name, comp_info in COMPARISONS.items():
            print(f"  {comp_name}:", end=" ")
            deg_df = run_deg_ttest(adata, comp_info)
            if deg_df is None:
                print("  insufficient cells")
                row[f"{comp_name}_rank"] = None; row[f"{comp_name}_nes"] = None; row[f"{comp_name}_inverted"] = None
                continue

            enr_df = run_gsea(deg_df)
            if enr_df is None:
                print("  GSEA failed")
                row[f"{comp_name}_rank"] = None; row[f"{comp_name}_nes"] = None; row[f"{comp_name}_inverted"] = None
                continue

            df_sorted = enr_df.sort_values("NES", ascending=False).reset_index(drop=True)
            nfkb_mask = df_sorted["Term"] == PATHWAY_NAME
            if nfkb_mask.sum() > 0:
                rank = df_sorted[nfkb_mask].index[0] + 1
                nes = df_sorted[nfkb_mask]["NES"].values[0]
                inverted = MAX_PATHWAYS - rank + 1
                print(f"rank={rank}/50, NES={nes:.4f} (inverted={inverted})")
            else:
                rank, nes, inverted = None, None, None
                print("NF-κB not found")

            row[f"{comp_name}_rank"] = rank
            row[f"{comp_name}_nes"] = nes
            row[f"{comp_name}_inverted"] = inverted
        results.append(row)

    df = pd.DataFrame(results)

    # Save
    out_csv = OUTPUT_DIR / "nfkb_rankings_13types_ttest_reproduce.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nSaved: {out_csv}")
    print(df.to_string(index=False))

    # Compare with previous backup (optional — only if available locally)
    backup_csv = GSEA_DIR / "nfkb_rankings_13types_mast_nr_vs_r_flipped.csv"
    if backup_csv.exists():
        old = pd.read_csv(backup_csv)
        print("\n\n=== COMPARISON WITH ORIGINAL ===")
        for _, r in df.iterrows():
            ct = r['cell_type']
            old_row = old[old['cell_type'] == ct]
            if len(old_row) == 0:
                continue
            o = old_row.iloc[0]
            pre_match = "✓" if r['pre_nes'] is not None and o['pre_nes'] is not None and \
                        abs(float(r['pre_nes']) - float(o['pre_nes'])) < 0.01 else "✗"
            post_match = "✓" if r['post_nes'] is not None and o['post_nes'] is not None and \
                         abs(float(r['post_nes']) - float(o['post_nes'])) < 0.01 else "✗"
            print(f"  {ct:20s}  pre: {str(r.get('pre_nes','N/A')):>8s} vs {o['pre_nes']:>8.4f} {pre_match}  "
                  f"post: {str(r.get('post_nes','N/A')):>8s} vs {o['post_nes']:>8.4f} {post_match}")


if __name__ == "__main__":
    main()
