#!/usr/bin/env python3
"""
Panel G: GSEA enrichment plots for NF-κB pathway — Fibroblast + Epithelial.
1×2 stacked layout. Direct subplot drawing (no intermediate PNG).
4× scaling method for Nature Cancer.

Note: Uses on-the-fly t-test DEGs (not MAST) with min_size=5.
This differs from the MAST-based GSEA used in Figs 5A/5F (min_size=15).
See Methods for details.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gseapy as gp
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *

BASE_DIR = Path(__file__).parent

# 4× scaling
DPI = 300
SCALE = 4
CM_TO_INCH = 1 / 2.54
PANEL_WIDTH_CM = 6.0 * SCALE
PANEL_HEIGHT_CM = 7.0 * SCALE

PATHWAY_NAME = "TNF-alpha Signaling via NF-kB"

CELL_TYPE_CONFIG = {
    "Fibroblast": {"path": FIBROBLAST_H5AD, "label": "Fibroblast"},
    "Epithelial": {"path": EPITHELIAL_H5AD, "label": "Epithelial"},
}

COMPARISON = {"column": "stomach_post_grouping", "group1": "Responsed", "group2": "No-response"}


def run_deg_ttest(adata, comparison_info, max_cells=5000):
    col = comparison_info["column"]
    g1, g2 = comparison_info["group1"], comparison_info["group2"]
    if col not in adata.obs.columns:
        return None
    valid_mask = adata.obs[col].isin([g1, g2])
    adata_sub = adata[valid_mask].copy()
    n_g1 = (adata_sub.obs[col] == g1).sum()
    n_g2 = (adata_sub.obs[col] == g2).sum()
    print(f"    R: {n_g1}, NR: {n_g2}")
    if n_g1 < 10 or n_g2 < 10:
        return None
    if adata_sub.n_obs > max_cells:
        sc.pp.subsample(adata_sub, n_obs=max_cells, random_state=42)
    try:
        sc.tl.rank_genes_groups(adata_sub, groupby=col, groups=[g2], reference=g1,
                                method='t-test_overestim_var', pts=True)
        return sc.get.rank_genes_groups_df(adata_sub, group=g2)
    except Exception as e:
        print(f"    DEG error: {e}")
        return None


def run_gsea(deg_df, cell_type):
    if deg_df is None or len(deg_df) == 0:
        return None
    deg_df = deg_df.copy()
    deg_df['rank_metric'] = deg_df['logfoldchanges'] * (-np.log10(deg_df['pvals'].clip(lower=1e-300)))
    rnk = deg_df.set_index('names')['rank_metric'].dropna()
    rnk = rnk[~rnk.index.duplicated(keep='first')].sort_values(ascending=False)
    if len(rnk) < 15:
        return None
    gsea_outdir = BASE_DIR / f"gsea_{cell_type}"
    gsea_outdir.mkdir(parents=True, exist_ok=True)
    try:
        pre_res = gp.prerank(rnk=rnk, gene_sets='MSigDB_Hallmark_2020',
                             outdir=str(gsea_outdir), min_size=5, max_size=500,
                             permutation_num=1000, seed=42, verbose=False)
        return pre_res
    except Exception as e:
        print(f"    GSEA error: {e}")
        return None


def draw_gsea_subplot(ax_top, ax_bot, pre_res, label):
    """Draw GSEA enrichment on two axes (ES curve + gene hits)."""
    term_data = pre_res.results[PATHWAY_NAME]
    RES = term_data['RES']
    hits = term_data['hits']
    nes = term_data['nes']
    pval = term_data['pval']
    fdr = term_data['fdr']

    x = np.arange(len(RES))
    color = '#C62828'

    # Top: running ES
    ax_top.plot(x, RES, color=color, linewidth=1.5)
    ax_top.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
    ax_top.fill_between(x, 0, RES, where=(np.array(RES) >= 0), color='#EF9A9A', alpha=0.3)
    ax_top.fill_between(x, 0, RES, where=(np.array(RES) < 0), color='#90CAF9', alpha=0.3)

    peak_idx = np.argmax(np.abs(RES))
    ax_top.axvline(x=peak_idx, color='red', linestyle=':', linewidth=0.8, alpha=0.7)

    ax_top.set_ylabel('Enrichment Score', fontsize=6 * SCALE)
    ax_top.set_xlim(0, len(RES))
    ax_top.tick_params(axis='both', labelsize=5 * SCALE, width=1.0, length=4)
    ax_top.set_title(label, fontsize=7 * SCALE)

    stats_text = f"NES = {nes:.2f}\np = {pval:.3f}\nFDR = {fdr:.3f}"
    ax_top.text(0.98, 0.95, stats_text, transform=ax_top.transAxes,
                fontsize=5 * SCALE, va='top', ha='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray',
                          linewidth=0.5))

    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)
    ax_top.spines['left'].set_linewidth(1.0)
    ax_top.spines['bottom'].set_linewidth(1.0)

    # Bottom: gene hits
    for hit in hits:
        ax_bot.axvline(x=hit, color='black', linewidth=0.5, alpha=0.5)
    ax_bot.set_xlim(0, len(RES))
    ax_bot.set_xlabel('Gene Rank', fontsize=6 * SCALE)
    ax_bot.set_ylabel('Hits', fontsize=6 * SCALE)
    ax_bot.set_yticks([])
    ax_bot.tick_params(axis='x', labelsize=5 * SCALE, width=1.0, length=4)
    ax_bot.spines['top'].set_visible(False)
    ax_bot.spines['right'].set_visible(False)
    ax_bot.spines['left'].set_linewidth(1.0)
    ax_bot.spines['bottom'].set_linewidth(1.0)

    print(f"    {label}: NES={nes:.2f}, p={pval:.3f}, FDR={fdr:.3f}")


def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    print("=" * 60)
    print("Panel G: GSEA Enrichment — Fibroblast + Epithelial")
    print("  (t-test DEGs, min_size=5)")
    print("=" * 60)

    # Run GSEA for each cell type
    gsea_results = {}
    for cell_type, config in CELL_TYPE_CONFIG.items():
        print(f"\n=== {cell_type} ===")
        if not config["path"].exists():
            print(f"  Not found: {config['path']}")
            continue
        adata = sc.read_h5ad(config["path"])
        adata.obs_names_make_unique()
        print(f"  {adata.n_obs} cells")

        deg_df = run_deg_ttest(adata, COMPARISON)
        if deg_df is None:
            continue
        pre_res = run_gsea(deg_df, cell_type)
        if pre_res is not None and PATHWAY_NAME in pre_res.results:
            gsea_results[cell_type] = pre_res

    if not gsea_results:
        print("No GSEA results! Exiting.")
        return

    # Draw directly in subplots — 2 cell types × (ES + hits) = 4 rows
    n_types = len(gsea_results)
    fig, axes = plt.subplots(n_types * 2, 1,
                              figsize=(PANEL_WIDTH_CM * CM_TO_INCH, PANEL_HEIGHT_CM * CM_TO_INCH),
                              gridspec_kw={'height_ratios': [3, 1] * n_types})

    for idx, (cell_type, pre_res) in enumerate(gsea_results.items()):
        ax_top = axes[idx * 2]
        ax_bot = axes[idx * 2 + 1]
        label = CELL_TYPE_CONFIG[cell_type]["label"]
        draw_gsea_subplot(ax_top, ax_bot, pre_res, label)

    plt.tight_layout()

    output = BASE_DIR / 'panel_g_gsea_4types.png'
    plt.savefig(output, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.savefig(output.with_suffix('.svg'), bbox_inches='tight', facecolor='white')
    plt.savefig(output.with_suffix('.pdf'), format='pdf', dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"\nSaved: {output}")


if __name__ == "__main__":
    main()
