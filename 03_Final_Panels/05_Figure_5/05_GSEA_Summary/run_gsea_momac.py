#!/usr/bin/env python3
"""Run GSEA for MoMac using the same t-test method as Panel I (05_G)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import MOMAC_H5AD

import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = Path(__file__).parent
COMPARISON = {"column": "stomach_post_grouping", "group1": "Responsed", "group2": "No-response"}


def main():
    print("=" * 60)
    print("GSEA for MoMac — t-test method")
    print("=" * 60)

    adata = sc.read_h5ad(MOMAC_H5AD)
    adata.obs_names_make_unique()
    adata = adata[adata.obs['Sample site'] == 'Stomach'].copy()
    adata = adata[adata.obs['Treatment phase'] == 'Post'].copy()
    print(f"  Post-treatment stomach MoMac: {adata.n_obs}")

    col = COMPARISON["column"]
    g1, g2 = COMPARISON["group1"], COMPARISON["group2"]
    valid_mask = adata.obs[col].isin([g1, g2])
    adata_sub = adata[valid_mask].copy()
    n_g1 = (adata_sub.obs[col] == g1).sum()
    n_g2 = (adata_sub.obs[col] == g2).sum()
    print(f"  R: {n_g1}, NR: {n_g2}")

    max_cells = 5000
    if adata_sub.n_obs > max_cells:
        sc.pp.subsample(adata_sub, n_obs=max_cells, random_state=42)

    sc.tl.rank_genes_groups(adata_sub, groupby=col, groups=[g2], reference=g1,
                            method='t-test_overestim_var', pts=True)
    deg_df = sc.get.rank_genes_groups_df(adata_sub, group=g2)

    # Ranking
    deg_df['rank_metric'] = deg_df['logfoldchanges'] * (-np.log10(deg_df['pvals'].clip(lower=1e-300)))
    rnk = deg_df.set_index('names')['rank_metric'].dropna()
    rnk = rnk[~rnk.index.duplicated(keep='first')].sort_values(ascending=False)

    gsea_outdir = BASE_DIR / "gsea_MoMac"
    gsea_outdir.mkdir(parents=True, exist_ok=True)

    pre_res = gp.prerank(rnk=rnk, gene_sets='MSigDB_Hallmark_2020',
                         outdir=str(gsea_outdir), min_size=5, max_size=500,
                         permutation_num=1000, seed=42, verbose=False)

    res_df = pre_res.res2d.copy()
    res_df = res_df.rename(columns={
        'Term': 'Pathway', 'NOM p-val': 'NOM_pval', 'FDR q-val': 'FDR_qval',
    })
    res_df['CellType'] = 'MoMac'
    gsea_data_dir = BASE_DIR / 'gsea_data'
    gsea_data_dir.mkdir(parents=True, exist_ok=True)
    res_df[['CellType', 'Pathway', 'NES', 'NOM_pval', 'FDR_qval']].to_csv(
        gsea_data_dir / 'gsea_momac.csv', index=False)

    # Print all Hallmark pathways sorted by NES
    res_df['NES'] = res_df['NES'].astype(float)
    res_df['NOM_pval'] = res_df['NOM_pval'].astype(float)
    res_df = res_df.sort_values('NES', ascending=False)

    print(f"\n  All {len(res_df)} Hallmark pathways (sorted by NES):")
    for _, r in res_df.iterrows():
        sig = '*' if r['NOM_pval'] < 0.05 else ''
        print(f"    NES={float(r['NES']):+.3f}  p={float(r['NOM_pval']):.4f} {sig:2s} {r['Pathway']}")


if __name__ == '__main__':
    main()
