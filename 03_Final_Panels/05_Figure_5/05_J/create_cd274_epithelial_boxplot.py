#!/usr/bin/env python3
"""Panel K: CD274 Post-R vs Post-NR in Epithelial cells"""
import warnings, numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import EPITHELIAL_H5AD

warnings.filterwarnings('ignore')

H5AD = EPITHELIAL_H5AD
OUT_DIR = Path(__file__).parent
SCALE = 4; CM = 1 / 2.54

BOX_COLORS = {'R': '#2166AC', 'NR': '#B2182B'}
MIN_CELLS = 20

def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'svg.fonttype': 'none',
        'pdf.fonttype': 42, 'ps.fonttype': 42,
    })

    adata = sc.read_h5ad(H5AD)
    adata = adata[adata.obs['Sample site'] == 'Stomach'].copy()
    adata = adata[adata.obs['Treatment phase'] == 'Post'].copy()
    adata.obs['response'] = adata.obs['stomach_post_grouping'].map({
        'Responsed': 'R', 'No-response': 'NR'}).astype(str)
    adata = adata[adata.obs['response'].isin(['R', 'NR'])].copy()

    gene = 'CD274'
    gene_names = list(adata.raw.var_names) if adata.raw else list(adata.var_names)
    if adata.raw:
        idx = list(adata.raw.var_names).index(gene)
        expr = adata.raw.X[:, idx]
    else:
        idx = list(adata.var_names).index(gene)
        expr = adata.X[:, idx]
    if hasattr(expr, 'toarray'):
        expr = expr.toarray().flatten()
    adata.obs['value'] = np.array(expr).flatten()

    # Sample-level aggregation
    df = adata.obs[['sample', 'response', 'value']].copy()
    df['sample'] = df['sample'].astype(str)
    counts = df.groupby('sample', observed=True).size()
    valid = counts[counts >= MIN_CELLS].index
    df = df[df['sample'].isin(valid)]
    sample_df = df.groupby(['sample', 'response'], observed=True)['value'].mean().reset_index()

    r_data = sample_df[sample_df['response'] == 'R']['value'].values
    nr_data = sample_df[sample_df['response'] == 'NR']['value'].values

    fig, ax = plt.subplots(figsize=(3.5 * SCALE * CM, 5 * SCALE * CM))
    bp = ax.boxplot([r_data, nr_data], positions=[0, 1], widths=0.5,
                   patch_artist=True, showfliers=False,
                   medianprops=dict(color='black', linewidth=1.5 * SCALE),
                   whiskerprops=dict(linewidth=1.0 * SCALE),
                   capprops=dict(linewidth=1.0 * SCALE),
                   boxprops=dict(linewidth=1.0 * SCALE))
    bp['boxes'][0].set_facecolor(BOX_COLORS['R']); bp['boxes'][0].set_alpha(0.6)
    bp['boxes'][1].set_facecolor(BOX_COLORS['NR']); bp['boxes'][1].set_alpha(0.6)

    for k, (data, color) in enumerate(zip([r_data, nr_data], [BOX_COLORS['R'], BOX_COLORS['NR']])):
        jitter = np.random.uniform(-0.08, 0.08, len(data))
        ax.scatter([k] * len(data) + jitter, data, c=color, s=20 * SCALE,
                  edgecolors='white', linewidths=0.3 * SCALE, alpha=0.85, zorder=3)

    _, pval = stats.mannwhitneyu(nr_data, r_data, alternative='greater')
    p_str = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'
    is_star = pval < 0.05
    y_max = max(np.max(r_data), np.max(nr_data))
    y_range = y_max - min(np.min(r_data), np.min(nr_data))
    bracket_y = y_max + 0.10 * y_range
    ax.plot([0, 0, 1, 1], [bracket_y - 0.02 * y_range, bracket_y,
            bracket_y, bracket_y - 0.02 * y_range], color='black', linewidth=0.8 * SCALE)
    ax.text(0.5, bracket_y + 0.02 * y_range, p_str, ha='center',
            fontsize=(7 if is_star else 5) * SCALE,
            fontweight='bold' if is_star else 'normal')

    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"R\n(n={len(r_data)})", f"NR\n(n={len(nr_data)})"],
                       fontsize=5 * SCALE)
    ax.set_ylabel('', fontsize=5 * SCALE)
    ax.set_title('PD-L1 (CD274)\nEpithelial', fontsize=6 * SCALE, fontweight='normal')
    ax.tick_params(labelsize=5 * SCALE, width=1.0 * SCALE, length=4 * SCALE)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for sp in ['left', 'bottom']:
        ax.spines[sp].set_linewidth(1.0 * SCALE)

    plt.tight_layout()
    stem = 'cd274_epithelial_boxplot'
    for ext in ['svg', 'png']:
        kw = {'bbox_inches': 'tight', 'facecolor': 'white'}
        if ext == 'png':
            kw['dpi'] = 300
        fig.savefig(OUT_DIR / f'{stem}.{ext}', **kw)
    print(f"R n={len(r_data)}, NR n={len(nr_data)}, P={pval:.4f}")
    print(f"Saved: {stem}")
    plt.close()


if __name__ == '__main__':
    main()
