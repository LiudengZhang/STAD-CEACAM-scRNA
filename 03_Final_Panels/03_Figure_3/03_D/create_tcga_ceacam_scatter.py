#!/usr/bin/env python3
"""
Figure 3 Row 2 (TCGA): CEACAM5/6 vs CD274 and Tex ssGSEA in TCGA-STAD.
Two 1×2 panels matching Row 1 layout:
  - tcga_ceacam_cd274.png: CEACAM5/6 vs CD274 (matches Panel A)
  - tcga_ceacam_tex.png:   CEACAM5/6 vs Tex ssGSEA (matches Panel B)

CD274: from BayesPrism deconvolved epithelial expression (log2+1)
Tex: proper ssGSEA (weighted KS) on bulk FPKM with 18-gene Tex signature
CEACAM: from BayesPrism deconvolved epithelial expression (log2+1)

Uses tumor-only 407 samples and gene mapping from STAR raw files.
"""

import pandas as pd
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import TCGA_BAYESPRISM_EPI, TCGA_BULK_TUMOR_ONLY, TCGA_RAW_DIR, TCGA_CLINICAL

# ==============================================================================
# Configuration
# ==============================================================================
SCALE = 4
DPI = 300
# 1×2 panels matching Panel A/B dimensions
PANEL_WIDTH_CM = 8.0 * SCALE
PANEL_HEIGHT_CM = 4.0 * SCALE
CM_TO_INCH = 1 / 2.54

OUTPUT_DIR = Path(__file__).parent

# Data paths (from paths.py)
TCGA_EPI_EXPR = TCGA_BAYESPRISM_EPI

# Tex signature (18 genes)
TEX_GENES = [
    'PDCD1', 'HAVCR2', 'LAG3', 'TIGIT', 'CTLA4', 'TOX', 'ENTPD1', 'CXCL13',
    'LAYN', 'CD38', 'BATF', 'IRF4', 'PRDM1', 'TOX2', 'ITGAE', 'NR4A1', 'NR4A2', 'NR4A3',
]

# Survival dot colors
COLOR_DEAD = '#c0392b'      # dark red, filled
COLOR_ALIVE = '#888888'     # gray edge, unfilled
COLOR_MISSING = '#cccccc'   # light gray edge, unfilled


# ==============================================================================
# ssGSEA scoring (proper weighted KS, from original validation session)
# ==============================================================================
def ssgsea_score(expr_df, gene_set):
    """Proper ssGSEA: weighted Kolmogorov-Smirnov enrichment. expr_df: samples x genes."""
    gene_set = [g for g in gene_set if g in expr_df.columns]
    print(f"  ssGSEA: {len(gene_set)}/{len(TEX_GENES)} Tex genes found")

    scores = []
    for sample in expr_df.index:
        profile = expr_df.loc[sample].dropna()
        ranked = profile.rank(ascending=True)
        sorted_idx = ranked.sort_values().index
        in_set_sorted = sorted_idx.isin(gene_set)
        ranks_sorted = ranked[sorted_idx]

        hits = in_set_sorted.astype(float)
        hw = hits * np.abs(ranks_sorted)
        hs = hw.sum()
        if hs == 0:
            scores.append(0)
            continue
        hr = np.cumsum(hw / hs)

        mw = (1 - hits)
        ms = mw.sum()
        if ms == 0:
            scores.append(0)
            continue
        mr = np.cumsum(mw / ms)
        scores.append((hr - mr).sum())

    return pd.Series(scores, index=expr_df.index)


# ==============================================================================
# Data loading
# ==============================================================================
def load_data():
    # 1. Deconvolved epithelial expression (samples × genes)
    print("[1/2] Loading TCGA BayesPrism epithelial expression...")
    epi = pd.read_csv(TCGA_EPI_EXPR, sep='\t', index_col=0)
    print(f"  Shape: {epi.shape}")

    # 2. Tumor-only FPKM (407 samples) + gene mapping from STAR raw files
    print("[2/2] Loading TCGA tumor-only FPKM for ssGSEA...")
    fpkm = pd.read_csv(TCGA_BULK_TUMOR_ONLY, sep='\t', index_col=0)
    fpkm = fpkm[fpkm.index.str.startswith('ENSG')]
    print(f"  FPKM: {fpkm.shape[0]:,} genes × {fpkm.shape[1]} samples")

    # Gene mapping from raw STAR annotation
    raw_file = os.listdir(TCGA_RAW_DIR)[0]
    ann = pd.read_csv(TCGA_RAW_DIR / raw_file, sep='\t', skiprows=1,
                       usecols=['gene_id', 'gene_name']).dropna()
    ens2sym = dict(zip(ann['gene_id'], ann['gene_name']))

    fpkm.index = fpkm.index.map(lambda x: ens2sym.get(x, x))
    fpkm = fpkm[~fpkm.index.duplicated(keep='first')]

    # Samples × genes, log2(FPKM+1)
    fpkm_log = np.log2(fpkm.T + 1)
    print(f"  After mapping: {fpkm_log.shape[1]:,} genes, {fpkm_log.shape[0]} samples")

    # ssGSEA
    print("  Running ssGSEA...")
    tex_scores = ssgsea_score(fpkm_log, TEX_GENES)

    # Align: epi (410 samples) ∩ FPKM (407 tumor-only)
    common = sorted(set(epi.index) & set(tex_scores.index))
    print(f"  Common samples: {len(common)}")

    df = pd.DataFrame({
        'sample': common,
        'CEACAM5': np.log2(epi.loc[common, 'CEACAM5'].values + 1),
        'CEACAM6': np.log2(epi.loc[common, 'CEACAM6'].values + 1),
        'CD274': np.log2(epi.loc[common, 'CD274'].values + 1),
        'tex_ssgsea_raw': tex_scores.loc[common].values,
    })

    # Z-score the Tex ssGSEA (z-score)s
    mu = df['tex_ssgsea_raw'].mean()
    sd = df['tex_ssgsea_raw'].std()
    df['tex_ssgsea'] = (df['tex_ssgsea_raw'] - mu) / sd
    print(f"  Tex ssGSEA z-scored: mean={mu:.2f}, sd={sd:.2f}")

    # Merge clinical survival data
    print("  Merging clinical vital_status...")
    clin = pd.read_csv(TCGA_CLINICAL, sep='\t', usecols=['submitter_id', 'vital_status'])
    clin = clin.drop_duplicates(subset='submitter_id')
    df = df.merge(clin, left_on='sample', right_on='submitter_id', how='left')
    df['vital_status'] = df['vital_status'].fillna('Missing')
    print(f"  Vital status counts: {df['vital_status'].value_counts().to_dict()}")

    return df


# ==============================================================================
# Plotting
# ==============================================================================
def plot_scatter(ax, x, y, vital, title, xlabel, ylabel, fontscale):
    r_val, p_val = stats.spearmanr(x, y)

    # Plot dots by vital_status: Dead filled red, Alive open gray, Missing open light gray
    for status, fc, ec, label in [
        ('Dead',    COLOR_DEAD,  COLOR_DEAD,    'Dead'),
        ('Alive',   'none',      COLOR_ALIVE,   'Alive'),
        ('Missing', 'none',      COLOR_MISSING, None),
    ]:
        mask = vital == status
        if mask.sum() == 0:
            continue
        ax.scatter(x[mask], y[mask], facecolors=fc, edgecolors=ec,
                   s=20*fontscale, alpha=0.85, linewidths=0.5*fontscale, zorder=3)

    # Regression line (above dots)
    valid = np.isfinite(x) & np.isfinite(y)
    xv, yv = np.array(x)[valid], np.array(y)[valid]
    if len(xv) >= 3 and np.std(xv) > 0:
        slope, intercept = np.polyfit(xv, yv, 1)
        x_line = np.linspace(xv.min(), xv.max(), 100)
        ax.plot(x_line, slope * x_line + intercept, 'k--', linewidth=0.8*fontscale, alpha=0.6, zorder=5)

    # Stats
    p_str = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'

    ax.set_title(f'{title}\nρ = {r_val:.2f}, {p_str}',
                 fontsize=6.5 * fontscale, fontweight='normal', linespacing=1.4)
    ax.set_xlabel(xlabel, fontsize=6 * fontscale)
    ax.set_ylabel(ylabel, fontsize=6 * fontscale)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for spine in ['bottom', 'left']:
        ax.spines[spine].set_linewidth(0.5*fontscale)
    ax.tick_params(axis='both', labelsize=5 * fontscale, width=0.5*fontscale, length=3*fontscale)
    ax.set_box_aspect(1)

    # Compact legend (Dead filled red, Alive open gray)
    h_dead = mlines.Line2D([], [], marker='o', color='none', markerfacecolor=COLOR_DEAD,
                           markeredgecolor=COLOR_DEAD, markersize=4*fontscale/SCALE,
                           label='Dead')
    h_alive = mlines.Line2D([], [], marker='o', color='none', markerfacecolor='none',
                            markeredgecolor=COLOR_ALIVE, markersize=4*fontscale/SCALE,
                            markeredgewidth=0.5, label='Alive')
    ax.legend(handles=[h_dead, h_alive], fontsize=4.5*fontscale, loc='upper left',
              frameon=False, handletextpad=0.3, borderpad=0.2)

    return r_val, p_val


def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 6 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    df = load_data()
    n = len(df)
    fig_w = PANEL_WIDTH_CM * CM_TO_INCH
    fig_h = PANEL_HEIGHT_CM * CM_TO_INCH

    # Panel D: CEACAM vs Tex ssGSEA only (filter zero-CEACAM samples per subplot)
    panel = {
        'name': 'tcga_ceacam_tex',
        'combos': [
            ('CEACAM5', 'tex_ssgsea', 'Epi. CEACAM5 (log$_2$+1)', 'Tex ssGSEA (z-score)'),
            ('CEACAM6', 'tex_ssgsea', 'Epi. CEACAM6 (log$_2$+1)', 'Tex ssGSEA (z-score)'),
        ],
    }

    results = []
    fig, axes = plt.subplots(1, 2, figsize=(fig_w, fig_h))
    for col_idx, (xcol, ycol, xlabel, ylabel) in enumerate(panel['combos']):
        ax = axes[col_idx]
        # Remove samples with zero CEACAM expression
        mask = df[xcol] > 0
        sub = df[mask]
        n_sub = len(sub)
        cohort_title = f'TCGA-STAD (n = {n_sub})'
        r_val, p_val = plot_scatter(ax, sub[xcol].values, sub[ycol].values,
                                    sub['vital_status'].values,
                                    cohort_title, xlabel, ylabel, SCALE)
        results.append({'x': xcol, 'y': ycol, 'r': r_val, 'p': p_val, 'n': n_sub})

    plt.tight_layout()

    out = OUTPUT_DIR / f"{panel['name']}.png"
    fig.savefig(out, dpi=DPI, bbox_inches='tight', facecolor='white')
    fig.savefig(str(out).replace('.png', '.svg'), format='svg', bbox_inches='tight', facecolor='white')
    print(f"Saved: {out}")
    plt.close()

    # Summary
    print(f"\n{'='*65}")
    print(f"TCGA-STAD CORRELATION SUMMARY (Spearman)")
    print(f"{'='*65}")
    print(f"{'X':<12} {'Y':<18} {'n':>5} {'r':>6} {'P':>12} {'Sig':>5}")
    print(f"{'-'*65}")
    for r in results:
        sig = '***' if r['p'] <= 0.001 else '**' if r['p'] <= 0.01 else '*' if r['p'] <= 0.05 else 'ns'
        print(f"{r['x']:<12} {r['y']:<18} {r['n']:>5} {r['r']:>6.3f} {r['p']:>12.6f} {sig:>5}")
    print(f"{'='*65}")


if __name__ == '__main__':
    main()
