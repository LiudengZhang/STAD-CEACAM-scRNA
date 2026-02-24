#!/usr/bin/env python3
"""
Run GSEA for 5 cell types using the same t-test method as Panel I (05_G).
- CD4+ T, CD8+ T, DC: compute fresh GSEA
- Epithelial, Fibroblast: read existing CSVs from 05_G
Combine all into gsea_combined_5types.csv
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import TCD4_H5AD, TCD8_H5AD, DC_CELLS_H5AD

import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = Path(__file__).parent
COMPARISON = {"column": "stomach_post_grouping", "group1": "Responsed", "group2": "No-response"}

# Existing GSEA CSVs
EXISTING_GSEA = {
    "Epithelial": Path(__file__).resolve().parents[1] / "05_G/gsea_Epithelial/gseapy.gene_set.prerank.report.csv",
    "Fibroblast": Path(__file__).resolve().parents[1] / "05_G/gsea_Fibroblast/gseapy.gene_set.prerank.report.csv",
}

# Cell types to compute fresh
COMPUTE_GSEA = {
    "CD4+ T": TCD4_H5AD,
    "CD8+ T": TCD8_H5AD,
    "DC": DC_CELLS_H5AD,
}


def run_deg_ttest(adata, max_cells=5000):
    col = COMPARISON["column"]
    g1, g2 = COMPARISON["group1"], COMPARISON["group2"]
    valid_mask = adata.obs[col].isin([g1, g2])
    adata_sub = adata[valid_mask].copy()
    n_g1 = (adata_sub.obs[col] == g1).sum()
    n_g2 = (adata_sub.obs[col] == g2).sum()
    print(f"    R: {n_g1}, NR: {n_g2}")
    if n_g1 < 10 or n_g2 < 10:
        print(f"    Too few cells, skipping")
        return None
    if adata_sub.n_obs > max_cells:
        sc.pp.subsample(adata_sub, n_obs=max_cells, random_state=42)
    sc.tl.rank_genes_groups(adata_sub, groupby=col, groups=[g2], reference=g1,
                            method='t-test_overestim_var', pts=True)
    return sc.get.rank_genes_groups_df(adata_sub, group=g2)


def run_gsea(deg_df, cell_type):
    deg_df = deg_df.copy()
    deg_df['rank_metric'] = deg_df['logfoldchanges'] * (-np.log10(deg_df['pvals'].clip(lower=1e-300)))
    rnk = deg_df.set_index('names')['rank_metric'].dropna()
    rnk = rnk[~rnk.index.duplicated(keep='first')].sort_values(ascending=False)
    if len(rnk) < 15:
        return None
    gsea_outdir = BASE_DIR / f"gsea_{cell_type.replace('+', 'plus').replace(' ', '_')}"
    gsea_outdir.mkdir(parents=True, exist_ok=True)
    pre_res = gp.prerank(rnk=rnk, gene_sets='MSigDB_Hallmark_2020',
                         outdir=str(gsea_outdir), min_size=5, max_size=500,
                         permutation_num=1000, seed=42, verbose=False)
    return pre_res


def main():
    print("=" * 60)
    print("GSEA for 5 cell types — t-test method (matching Panel I)")
    print("=" * 60)

    all_results = []

    # 1. Read existing CSVs for Epithelial and Fibroblast
    for cell_type, csv_path in EXISTING_GSEA.items():
        print(f"\n=== {cell_type} (existing CSV) ===")
        df = pd.read_csv(csv_path)
        df = df.rename(columns={
            'Term': 'Pathway',
            'NOM p-val': 'NOM_pval',
            'FDR q-val': 'FDR_qval',
        })
        df['CellType'] = cell_type
        all_results.append(df[['CellType', 'Pathway', 'NES', 'NOM_pval', 'FDR_qval']])
        print(f"    {len(df)} pathways loaded")

    # 2. Compute fresh GSEA for CD4+ T, CD8+ T, DC
    for cell_type, h5ad_path in COMPUTE_GSEA.items():
        print(f"\n=== {cell_type} (computing fresh) ===")
        adata = sc.read_h5ad(h5ad_path)
        adata.obs_names_make_unique()

        # Filter to stomach post-treatment
        adata = adata[adata.obs['Sample site'] == 'Stomach'].copy()
        adata = adata[adata.obs['Treatment phase'] == 'Post'].copy()
        print(f"    {adata.n_obs} cells")

        deg_df = run_deg_ttest(adata)
        if deg_df is None:
            print(f"    Skipping {cell_type}")
            continue

        pre_res = run_gsea(deg_df, cell_type)
        if pre_res is None:
            print(f"    GSEA failed for {cell_type}")
            continue

        # Extract results
        res_df = pre_res.res2d.copy()
        res_df = res_df.rename(columns={
            'Term': 'Pathway',
            'NOM p-val': 'NOM_pval',
            'FDR q-val': 'FDR_qval',
        })
        res_df['CellType'] = cell_type
        all_results.append(res_df[['CellType', 'Pathway', 'NES', 'NOM_pval', 'FDR_qval']])
        print(f"    {len(res_df)} pathways")

        # Print key pathways
        key_pathways = [
            'TNF-alpha Signaling via NF-kB', 'Inflammatory Response',
            'IL-6/JAK/STAT3 Signaling', 'Interferon Gamma Response',
            'Epithelial Mesenchymal Transition', 'Angiogenesis',
            'Hypoxia', 'Apoptosis'
        ]
        for pw in key_pathways:
            row = res_df[res_df['Pathway'] == pw]
            if len(row) > 0:
                r = row.iloc[0]
                print(f"    {pw}: NES={float(r['NES']):.3f}, p={float(r['NOM_pval']):.4f}")

    # 3. Combine all
    combined = pd.concat(all_results, ignore_index=True)
    combined['NES'] = combined['NES'].astype(float)
    combined['NOM_pval'] = combined['NOM_pval'].astype(float)
    combined['FDR_qval'] = combined['FDR_qval'].astype(float)

    out_csv = BASE_DIR / 'gsea_data' / 'gsea_combined_5types.csv'
    combined.to_csv(out_csv, index=False)
    print(f"\nSaved combined: {out_csv}")
    print(f"  Total rows: {len(combined)}")
    print(f"  Cell types: {combined['CellType'].unique().tolist()}")

    # Summary table for the 8 key pathways
    print("\n" + "=" * 80)
    print("Summary: 8 Key Pathways × 5 Cell Types")
    print("=" * 80)
    key_pathways = [
        'TNF-alpha Signaling via NF-kB', 'Inflammatory Response',
        'IL-6/JAK/STAT3 Signaling', 'Interferon Gamma Response',
        'Epithelial Mesenchymal Transition', 'Angiogenesis',
        'Hypoxia', 'Apoptosis'
    ]
    for pw in key_pathways:
        print(f"\n  {pw}:")
        for ct in ['Epithelial', 'Fibroblast', 'CD8+ T', 'CD4+ T', 'DC']:
            row = combined[(combined['CellType'] == ct) & (combined['Pathway'] == pw)]
            if len(row) > 0:
                r = row.iloc[0]
                sig = '*' if float(r['NOM_pval']) < 0.05 else ''
                print(f"    {ct:15s}: NES={float(r['NES']):+.3f}  p={float(r['NOM_pval']):.4f} {sig}")
            else:
                print(f"    {ct:15s}: not found")


if __name__ == '__main__':
    main()
