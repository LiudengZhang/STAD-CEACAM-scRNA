#!/usr/bin/env python3
"""
S6: CD274 (PD-L1) Expression — Remaining Cell Types
Post-R vs Post-NR boxplots for 9 cell types not shown in Figure 5 (DC in main figure).
Colors: Post-R #ffcfd2, Post-NR #f1c0e8
Mann-Whitney U test, sample-level aggregation (min 20 cells/sample).
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

SCALE = 4
DPI = 300
CM = 1 / 2.54
MIN_CELLS = 20

COLOR_R = '#ffcfd2'   # Post-R
COLOR_NR = '#f1c0e8'  # Post-NR

# S6_A through S6_J: cell type → (letter, h5ad_path, major_cell_type_name, display_name)
PANELS = [
    ('A', TCD8_H5AD,       'CD8+ T cells',             'CD8+ T cells'),
    ('B', TCD4_H5AD,       'CD4+ T cells',             'CD4+ T cells'),
    ('C', NK_CELLS_H5AD,   'NK cells',                 'NK cells'),
    ('D', B_CELLS_H5AD,    'B cells',                  'B cells'),
    ('E', None,             'Plasma cells',             'Plasma cells'),
    ('F', None,             'Mast cells',               'Mast cells'),
    ('G', NEUTROPHILS_H5AD, 'Neutrophils',             'Neutrophils'),
    ('H', ENDOTHELIAL_H5AD, 'Endothelial cells',       'Endothelial cells'),
    ('I', PERICYTE_H5AD,   'Pericyte',                 'Pericytes'),
]

OUT_DIR = Path(__file__).parent


def create_cd274_panel(letter, h5ad_path, cell_type_name, display_name, full_adata=None):
    """Create a single CD274 boxplot panel."""
    panel_dir = OUT_DIR / f'S6_{letter}'
    panel_dir.mkdir(exist_ok=True)

    print(f"\n  Panel S6_{letter}: {display_name}")

    # Load data
    if h5ad_path is not None and h5ad_path.exists():
        adata = sc.read_h5ad(h5ad_path)
    elif full_adata is not None:
        adata = full_adata[full_adata.obs['major_cell_type'] == cell_type_name].copy()
    else:
        print(f"    WARNING: No h5ad for {cell_type_name}, loading full dataset...")
        adata = sc.read_h5ad(FULL_DATASET_H5AD)
        adata = adata[adata.obs['major_cell_type'] == cell_type_name].copy()

    print(f"    Total cells: {adata.n_obs}")

    # Filter to stomach, post-treatment
    adata = adata[adata.obs['Sample site'] == 'Stomach'].copy()
    adata = adata[adata.obs['Treatment phase'] == 'Post'].copy()
    adata.obs['response'] = adata.obs['stomach_post_grouping'].map({
        'Responsed': 'R', 'No-response': 'NR'
    }).astype(str)
    adata = adata[adata.obs['response'].isin(['R', 'NR'])].copy()
    print(f"    Post-treatment stomach cells: {adata.n_obs}")

    if adata.n_obs == 0:
        print(f"    No cells available. Generating empty panel.")
        _save_empty_panel(panel_dir, letter, display_name)
        return

    # Get CD274 expression
    gene = 'CD274'
    if adata.raw is not None and gene in adata.raw.var_names:
        idx = list(adata.raw.var_names).index(gene)
        expr = adata.raw.X[:, idx]
    elif gene in adata.var_names:
        idx = list(adata.var_names).index(gene)
        expr = adata.X[:, idx]
    else:
        print(f"    CD274 not found in var_names!")
        _save_empty_panel(panel_dir, letter, display_name)
        return

    if hasattr(expr, 'toarray'):
        expr = expr.toarray().flatten()
    adata.obs['cd274_expr'] = np.array(expr).flatten()

    # Sample-level aggregation
    df = adata.obs[['sample', 'response', 'cd274_expr']].copy()
    df['sample'] = df['sample'].astype(str)
    counts = df.groupby('sample', observed=True).size()
    valid = counts[counts >= MIN_CELLS].index
    df = df[df['sample'].isin(valid)]
    sample_df = df.groupby(['sample', 'response'], observed=True)['cd274_expr'].mean().reset_index()

    r_data = sample_df[sample_df['response'] == 'R']['cd274_expr'].values
    nr_data = sample_df[sample_df['response'] == 'NR']['cd274_expr'].values
    print(f"    R samples: n={len(r_data)}, NR samples: n={len(nr_data)}")

    # Mann-Whitney
    if len(r_data) >= 2 and len(nr_data) >= 2:
        _, pval = mannwhitneyu(nr_data, r_data, alternative='two-sided')
    else:
        pval = np.nan
    print(f"    P-value: {pval:.4f}" if not np.isnan(pval) else "    P-value: N/A")

    # Plot
    fig, ax = plt.subplots(figsize=(3.5 * SCALE * CM, 5 * SCALE * CM))

    bp = ax.boxplot([r_data, nr_data], positions=[0, 1], widths=0.5,
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color='black', linewidth=1.5 * SCALE),
                    whiskerprops=dict(linewidth=1.0 * SCALE),
                    capprops=dict(linewidth=1.0 * SCALE),
                    boxprops=dict(linewidth=1.0 * SCALE))
    bp['boxes'][0].set_facecolor(COLOR_R)
    bp['boxes'][0].set_alpha(0.7)
    bp['boxes'][1].set_facecolor(COLOR_NR)
    bp['boxes'][1].set_alpha(0.7)

    # Jitter points
    rng = np.random.default_rng(42)
    for k, (data, color) in enumerate(zip([r_data, nr_data], [COLOR_R, COLOR_NR])):
        if len(data) > 0:
            jitter = rng.uniform(-0.08, 0.08, len(data))
            ax.scatter([k] * len(data) + jitter, data, c=color, s=20 * SCALE,
                       edgecolors='black', linewidths=0.3 * SCALE, alpha=0.85, zorder=3)

    # Significance bracket
    all_vals = np.concatenate([r_data, nr_data]) if (len(r_data) > 0 and len(nr_data) > 0) else np.array([0])
    y_max = np.max(all_vals)
    y_range = np.max(all_vals) - np.min(all_vals) if len(all_vals) > 1 else 0.1
    if y_range == 0:
        y_range = 0.1
    bh = y_max + 0.10 * y_range
    ax.plot([0, 0, 1, 1], [bh - 0.02 * y_range, bh, bh, bh - 0.02 * y_range],
            color='black', linewidth=0.8 * SCALE)

    if np.isnan(pval):
        p_str = 'N/A'
    else:
        p_str = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'
    ax.text(0.5, bh + 0.02 * y_range, p_str, ha='center',
            fontsize=5 * SCALE)

    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"Post-R\n(n={len(r_data)})", f"Post-NR\n(n={len(nr_data)})"],
                        fontsize=5 * SCALE)
    ax.set_ylabel('PD-L1 (CD274)\nExpression', fontsize=5 * SCALE)
    ax.set_title(f'PD-L1 (CD274)\n{display_name}', fontsize=6 * SCALE, fontweight='normal')

    ax.tick_params(labelsize=5 * SCALE, width=1.0 * SCALE, length=4 * SCALE)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for sp in ['left', 'bottom']:
        ax.spines[sp].set_linewidth(1.0 * SCALE)

    ax.set_ylim(ax.get_ylim()[0], bh + 0.15 * y_range)

    plt.tight_layout()

    stem = f'panel_S6_{letter}'
    for ext in ['svg', 'png']:
        kw = {'bbox_inches': 'tight', 'facecolor': 'white'}
        if ext == 'png':
            kw['dpi'] = DPI
        fig.savefig(panel_dir / f'{stem}.{ext}', **kw)
    print(f"    Saved: {stem}")
    plt.close()


def _save_empty_panel(panel_dir, letter, display_name):
    """Save an empty panel with 'No data' text."""
    fig, ax = plt.subplots(figsize=(3.5 * SCALE * CM, 5 * SCALE * CM))
    ax.text(0.5, 0.5, f'PD-L1 (CD274)\n{display_name}\n\nInsufficient data',
            ha='center', va='center', fontsize=5 * SCALE,
            transform=ax.transAxes)
    ax.set_axis_off()
    plt.tight_layout()

    stem = f'panel_S6_{letter}'
    for ext in ['svg', 'png']:
        kw = {'bbox_inches': 'tight', 'facecolor': 'white'}
        if ext == 'png':
            kw['dpi'] = DPI
        fig.savefig(panel_dir / f'{stem}.{ext}', **kw)
    plt.close()


def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    print("=" * 60)
    print("S6: CD274 Remaining Cell Types")
    print("=" * 60)

    # Load full dataset once for cell types without dedicated h5ad
    print("Loading full dataset for Plasma/Mast cells...")
    full_adata = sc.read_h5ad(FULL_DATASET_H5AD)

    for letter, h5ad_path, cell_type_name, display_name in PANELS:
        create_cd274_panel(letter, h5ad_path, cell_type_name, display_name,
                          full_adata=full_adata)

    del full_adata
    import gc
    gc.collect()

    print("\nAll 9 CD274 panels complete!")


if __name__ == '__main__':
    main()
