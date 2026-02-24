#!/usr/bin/env python3
"""
Panel A: CEACAM5/6 vs CD274 correlation (1x2)
Primary cohort only (scRNA, ALL 32 stomach samples, 4-group coloring)
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *

# 4x scaling
SCALE = 4
DPI = 300
PANEL_WIDTH_CM = 8.0 * SCALE   # wide enough for 1x2
PANEL_HEIGHT_CM = 4.0 * SCALE  # single row
CM_TO_INCH = 1 / 2.54

EPI_RAW = EPITHELIAL_RAW_COUNTS_H5AD
OUTPUT_DIR = Path(__file__).parent

# Standard 4-group palette (blue=R, warm=NR; lighter=Pre, darker=Post)
COLOR_MAP_4GROUP = {
    'Pre-R':  '#74b9ff',
    'Pre-NR': '#e17055',
    'Post-R': '#0984e3',
    'Post-NR':'#d63031',
    'Other':  '#999999',
}



def load_primary():
    """Load primary cohort: ALL 32 stomach samples (pre + post)."""
    print("[Primary] Loading epithelial raw counts...")
    adata = sc.read_h5ad(EPI_RAW)

    # Filter: stomach only (all treatment phases)
    stomach_mask = adata.obs['Sample site'].astype(str).str.lower().str.contains('stomach')
    adata = adata[stomach_mask].copy()

    # Normalize raw counts
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Extract expression
    genes = ['CD274', 'CEACAM5', 'CEACAM6']
    expr = {}
    for g in genes:
        x = adata[:, g].X
        expr[g] = x.toarray().flatten() if hasattr(x, 'toarray') else x.flatten()

    df = pd.DataFrame(expr)
    df['sample'] = adata.obs['sample'].values
    df['treatment_phase'] = adata.obs['Treatment phase'].values
    df['pre_group'] = adata.obs['stomach_pre_grouping'].values
    df['post_group'] = adata.obs['stomach_post_grouping'].values

    # Sample-level aggregation
    agg = df.groupby('sample').agg({
        'CD274': 'mean', 'CEACAM5': 'mean', 'CEACAM6': 'mean',
        'treatment_phase': 'first',
        'pre_group': 'first',
        'post_group': 'first',
    }).reset_index()

    # 4-group assignment
    def assign_4group(row):
        phase = str(row['treatment_phase'])
        if phase == 'Pre':
            grp = str(row['pre_group'])
            if grp == 'Responsed': return 'Pre-R'
            if grp == 'No-response': return 'Pre-NR'
        elif phase == 'Post':
            grp = str(row['post_group'])
            if grp == 'Responsed': return 'Post-R'
            if grp == 'No-response': return 'Post-NR'
        return 'Other'

    agg['group'] = agg.apply(assign_4group, axis=1)

    # Exclude GC_1228: statistical outlier on CD274 (Grubbs' test G=5.43, P<0.05)
    agg = agg[agg['sample'] != 'GC_1228'].reset_index(drop=True)

    counts = agg['group'].value_counts()
    print(f"  {len(agg)} stomach samples (after outlier removal): {counts.to_dict()}")
    return agg


def load_tiger():
    """Load TIGER PRJEB25780: 78 ICB-treated GC (all pre-treatment)."""
    print("[TIGER] Loading expression...")
    expr = pd.read_csv(TIGER_EXPR, sep='\t', index_col=0)
    genes = ['CD274', 'CEACAM5', 'CEACAM6']

    tiger = expr.loc[genes].T.copy()
    tiger.columns = genes
    tiger = np.log2(tiger + 1)

    # Metadata
    meta = pd.read_csv(TIGER_META, sep='\t')
    id_col = None
    for col in meta.columns:
        if meta[col].dtype == object and meta[col].isin(tiger.index).sum() > 10:
            id_col = col
            break
    if id_col:
        meta_indexed = meta.set_index(id_col)
        tiger = tiger.join(meta_indexed['response'], how='left')

    def map_resp(val):
        val = str(val)
        if val in ['CR', 'PR']: return 'Pre-R'
        if val in ['SD', 'PD']: return 'Pre-NR'
        return 'Other'
    tiger['group'] = tiger['response'].map(map_resp)
    tiger = tiger.reset_index().rename(columns={'index': 'sample'})

    print(f"  {len(tiger)} samples — Pre-R:{(tiger.group=='Pre-R').sum()}, Pre-NR:{(tiger.group=='Pre-NR').sum()}")
    return tiger


def plot_scatter(ax, x, y, groups, color_map, draw_order, title, xlabel, ylabel, fontscale):
    """Generic scatter with Spearman stats."""
    r_val, p_val = stats.spearmanr(x, y)

    for grp in draw_order:
        mask = groups == grp
        if mask.sum() > 0:
            ax.scatter(x[mask], y[mask], c=color_map[grp], s=30*SCALE,
                       alpha=0.85, edgecolors='white', linewidths=0.3*SCALE,
                       label=grp, zorder=3)

    # Regression line
    valid = np.isfinite(x) & np.isfinite(y)
    xv, yv = x[valid], y[valid]
    if len(xv) >= 3 and xv.std() > 0:
        slope, intercept = np.polyfit(xv, yv, 1)
        x_line = np.linspace(xv.min(), xv.max(), 100)
        ax.plot(x_line, slope * x_line + intercept, 'k--', linewidth=0.8*SCALE, alpha=0.6, zorder=2)

    # Stats text — asterisk notation
    p_str = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'

    ax.set_title(f'{title}\nρ = {r_val:.2f}, {p_str}',
                 fontsize=6.5 * fontscale, fontweight='normal', linespacing=1.4)

    ax.set_xlabel(xlabel, fontsize=6 * fontscale)
    ax.set_ylabel(ylabel, fontsize=6 * fontscale)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for spine in ['bottom', 'left']:
        ax.spines[spine].set_linewidth(0.5*SCALE)
    ax.tick_params(axis='both', labelsize=5 * fontscale, width=0.5*SCALE, length=3*SCALE)

    ax.set_box_aspect(1)


def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 6 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    primary = load_primary()

    # 1x2 layout (primary cohort only)
    fig_w = PANEL_WIDTH_CM * CM_TO_INCH
    fig_h = PANEL_HEIGHT_CM * CM_TO_INCH
    fig, axes = plt.subplots(1, 2, figsize=(fig_w, fig_h))

    # Draw order: Other first (background), then colored groups on top
    primary_order = ['Other', 'Pre-R', 'Pre-NR', 'Post-R', 'Post-NR']

    for col_idx, ceacam in enumerate(['CEACAM5', 'CEACAM6']):
        plot_scatter(
            ax=axes[col_idx],
            x=primary[ceacam].values.astype(float),
            y=primary['CD274'].values.astype(float),
            groups=primary['group'].values,
            color_map=COLOR_MAP_4GROUP,
            draw_order=primary_order,
            title='Primary Cohort (scRNA-seq)',
            xlabel=f'{ceacam} (mean log expr.)',
            ylabel='PD-L1 (CD274) (mean log expr.)',
            fontscale=SCALE,
        )
        axes[col_idx].set_ylim(0, 0.03)

    plt.tight_layout()

    out = OUTPUT_DIR / 'ceacam_cd274_scatter.png'
    fig.savefig(out, dpi=DPI, bbox_inches='tight', facecolor='white')
    fig.savefig(str(out).replace('.png', '.pdf'), dpi=DPI, bbox_inches='tight', facecolor='white')
    fig.savefig(str(out).replace('.png', '.svg'), format='svg', bbox_inches='tight', facecolor='white')
    print(f"\nSaved: {out}")
    plt.close()


if __name__ == '__main__':
    main()
