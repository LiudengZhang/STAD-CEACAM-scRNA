#!/usr/bin/env python3
"""
S2_I: Liver C2_CEACAM6 Epithelial Cluster Proportion — R vs NR
Style matches Figure 2 02_E (no jitter, thin lines).
Uses liver_pre_grouping column for R/NR labels.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import EPITHELIAL_H5AD

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
CM_TO_INCH = 1 / 2.54
PANEL_WIDTH_CM = 3.2 * 4
PANEL_HEIGHT_CM = 3.0 * 4

# Colors matching Fig 2
COLORS = {'R': '#0072B2', 'NR': '#D55E00'}
MEDIAN_COLORS = {'R': '#005689', 'NR': '#A34700'}

OUT_DIR = Path(__file__).parent


def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 6 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    print("S2_I: Liver C2_CEACAM6 Proportion R vs NR")
    print("=" * 60)

    # Load epithelial h5ad
    print("Loading epithelial h5ad...")
    adata = sc.read_h5ad(EPITHELIAL_H5AD)
    print(f"  Total epithelial cells: {adata.n_obs}")

    # Filter to liver tissue
    adata = adata[adata.obs['Sample site'] == 'Liver'].copy()
    print(f"  Liver epithelial cells: {adata.n_obs}")

    # Map R/NR using liver_pre_grouping
    adata.obs['response'] = adata.obs['liver_pre_grouping'].map({
        'Responsed': 'R', 'No-response': 'NR'
    }).astype(str)
    adata = adata[adata.obs['response'].isin(['R', 'NR'])].copy()
    print(f"  After R/NR filtering: {adata.n_obs}")

    # Identify C2_CEACAM6 cells
    c2_mask = adata.obs['minor_cell_state'].str.contains('C2_', na=False) & \
              adata.obs['minor_cell_state'].str.contains('CEACAM', na=False)
    if c2_mask.sum() == 0:
        c2_mask = adata.obs['minor_cell_state'].str.contains('CEACAM', na=False)
    print(f"  C2_CEACAM cells: {c2_mask.sum()}")

    # Sample-level proportion
    adata.obs['is_c2'] = c2_mask.astype(int)
    sample_df = adata.obs.groupby(['sample', 'response'], observed=True).agg(
        n_total=('is_c2', 'size'),
        n_c2=('is_c2', 'sum')
    ).reset_index()
    sample_df['proportion'] = (sample_df['n_c2'] / sample_df['n_total']) * 100

    r_data = sample_df[sample_df['response'] == 'R']['proportion'].values
    nr_data = sample_df[sample_df['response'] == 'NR']['proportion'].values
    print(f"  R samples: n={len(r_data)}, NR samples: n={len(nr_data)}")

    # Mann-Whitney U test
    if len(r_data) >= 2 and len(nr_data) >= 2:
        _, pval = mannwhitneyu(r_data, nr_data, alternative='two-sided')
    else:
        pval = np.nan
    print(f"  P-value: {pval:.4f}" if not np.isnan(pval) else "  P-value: N/A")

    # Plot — Fig 2 style (no jitter, thin lines)
    fig, ax = plt.subplots(figsize=(PANEL_WIDTH_CM * CM_TO_INCH, PANEL_HEIGHT_CM * CM_TO_INCH))

    data = [r_data, nr_data]
    bp = ax.boxplot(data, positions=[1, 2], widths=0.6, patch_artist=True,
                    boxprops=dict(linewidth=0.5 * SCALE),
                    whiskerprops=dict(color='black', linewidth=0.5 * SCALE),
                    capprops=dict(color='black', linewidth=0.5 * SCALE),
                    flierprops=dict(marker='o', markerfacecolor='white', markersize=4,
                                   markeredgecolor='black', markeredgewidth=0.5 * SCALE),
                    medianprops=dict(linewidth=0.8 * SCALE))

    bp['boxes'][0].set_facecolor(COLORS['R'])
    bp['boxes'][1].set_facecolor(COLORS['NR'])
    bp['medians'][0].set_color(MEDIAN_COLORS['R'])
    bp['medians'][1].set_color(MEDIAN_COLORS['NR'])

    # Significance bracket
    all_vals = np.concatenate([r_data, nr_data])
    y_max = np.max(all_vals) if len(all_vals) > 0 else 1
    y_bracket = y_max * 1.15

    ax.plot([1, 1, 2, 2], [y_bracket, y_bracket * 1.05, y_bracket * 1.05, y_bracket],
            'k-', linewidth=0.5 * SCALE)

    if np.isnan(pval):
        p_str = 'N/A'
    else:
        p_str = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'
    ax.text(1.5, y_bracket * 1.08, p_str, ha='center', va='bottom', fontsize=6 * SCALE)

    ax.set_title('CEACAM5/6\nEpithelial', fontsize=6.5 * SCALE, fontweight='normal')
    ax.set_ylabel('Proportion (%)', fontsize=6 * SCALE)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['R', 'NR'], fontsize=6 * SCALE)
    ax.set_ylim(0, y_max * 1.40)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.5 * SCALE)
    ax.spines['bottom'].set_linewidth(0.5 * SCALE)
    ax.tick_params(axis='both', labelsize=5 * SCALE, width=0.5 * SCALE, length=3 * SCALE)
    fig.subplots_adjust(left=0.28, right=0.95, top=0.85, bottom=0.15)

    stem = 'panel_S2_I'
    for ext in ['svg', 'png']:
        kw = {'bbox_inches': 'tight', 'facecolor': 'white'}
        if ext == 'png':
            kw['dpi'] = DPI
        fig.savefig(OUT_DIR / f'{stem}.{ext}', **kw)
    print(f"Saved: {stem}")
    plt.close()


if __name__ == '__main__':
    main()
