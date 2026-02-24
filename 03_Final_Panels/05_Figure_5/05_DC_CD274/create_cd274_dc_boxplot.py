#!/usr/bin/env python3
"""Panel M: PD-L1 (CD274) boxplot — DC cells, Post-R vs Post-NR."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import DC_CELLS_H5AD

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

COLOR_R = '#2166AC'
COLOR_NR = '#B2182B'

OUT_DIR = Path(__file__).parent


def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    print("=" * 60)
    print("Panel M: CD274 DC cells — Post-R vs Post-NR")
    print("=" * 60)

    adata = sc.read_h5ad(DC_CELLS_H5AD)
    print(f"  Total cells: {adata.n_obs}")

    # Filter to stomach post-treatment
    adata = adata[adata.obs['Sample site'] == 'Stomach'].copy()
    adata = adata[adata.obs['Treatment phase'] == 'Post'].copy()
    adata.obs['response'] = adata.obs['stomach_post_grouping'].map({
        'Responsed': 'R', 'No-response': 'NR'
    }).astype(str)
    adata = adata[adata.obs['response'].isin(['R', 'NR'])].copy()
    print(f"  Post-treatment stomach cells: {adata.n_obs}")

    # Get CD274 expression
    gene = 'CD274'
    if adata.raw is not None and gene in adata.raw.var_names:
        idx = list(adata.raw.var_names).index(gene)
        expr = adata.raw.X[:, idx]
    elif gene in adata.var_names:
        idx = list(adata.var_names).index(gene)
        expr = adata.X[:, idx]
    else:
        print(f"  CD274 not found!")
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
    print(f"  R samples: n={len(r_data)}, NR samples: n={len(nr_data)}")

    # Mann-Whitney (one-sided: NR > R)
    if len(r_data) >= 2 and len(nr_data) >= 2:
        _, pval = mannwhitneyu(nr_data, r_data, alternative='greater')
    else:
        pval = np.nan
    print(f"  P-value (one-sided NR>R): {pval:.4f}" if not np.isnan(pval) else "  P-value: N/A")

    # Plot
    fig, ax = plt.subplots(figsize=(3.5 * SCALE * CM, 5 * SCALE * CM))

    bp = ax.boxplot([r_data, nr_data], positions=[0, 1], widths=0.5,
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color='black', linewidth=1.5 * SCALE),
                    whiskerprops=dict(linewidth=1.0 * SCALE),
                    capprops=dict(linewidth=1.0 * SCALE),
                    boxprops=dict(linewidth=1.0 * SCALE))
    bp['boxes'][0].set_facecolor(COLOR_R)
    bp['boxes'][0].set_alpha(0.6)
    bp['boxes'][1].set_facecolor(COLOR_NR)
    bp['boxes'][1].set_alpha(0.6)

    # Jitter points
    rng = np.random.default_rng(42)
    for k, (data, color) in enumerate(zip([r_data, nr_data], [COLOR_R, COLOR_NR])):
        if len(data) > 0:
            jitter = rng.uniform(-0.08, 0.08, len(data))
            ax.scatter([k] * len(data) + jitter, data, c=color, s=20 * SCALE,
                       edgecolors='white', linewidths=0.3 * SCALE, alpha=0.85, zorder=3)

    # Significance bracket
    all_vals = np.concatenate([r_data, nr_data])
    y_max = np.max(all_vals)
    y_range = np.max(all_vals) - np.min(all_vals)
    if y_range == 0:
        y_range = 0.1
    bh = y_max + 0.10 * y_range
    ax.plot([0, 0, 1, 1], [bh - 0.02 * y_range, bh, bh, bh - 0.02 * y_range],
            color='black', linewidth=0.8 * SCALE)

    p_str = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'
    is_star = pval < 0.05
    ax.text(0.5, bh + 0.02 * y_range, p_str, ha='center',
            fontsize=(7 if is_star else 5) * SCALE,
            fontweight='bold' if is_star else 'normal')

    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"R\n(n={len(r_data)})", f"NR\n(n={len(nr_data)})"],
                        fontsize=5 * SCALE)
    ax.set_ylabel('PD-L1 (CD274)\nExpression', fontsize=5 * SCALE)
    ax.set_title('PD-L1 (CD274)\nDC cells', fontsize=6 * SCALE, fontweight='normal')

    ax.tick_params(labelsize=5 * SCALE, width=1.0 * SCALE, length=4 * SCALE)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for sp in ['left', 'bottom']:
        ax.spines[sp].set_linewidth(1.0 * SCALE)

    ax.set_ylim(ax.get_ylim()[0], bh + 0.15 * y_range)

    plt.tight_layout()

    stem = 'cd274_dc_boxplot'
    for ext in ['svg', 'png']:
        kw = {'bbox_inches': 'tight', 'facecolor': 'white'}
        if ext == 'png':
            kw['dpi'] = DPI
        fig.savefig(OUT_DIR / f'{stem}.{ext}', **kw)
    print(f"  Saved: {stem}")
    plt.close()


if __name__ == '__main__':
    main()
