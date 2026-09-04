#!/usr/bin/env python3
"""
Figure 3 panel E, RESTYLED (Version B) - CEACAM5/6 vs Tex ssGSEA in TCGA-STAD.

NOTE THE PANEL LETTER, AND THE FILENAME. This directory is `03_D` but
PROVENANCE.csv says it holds printed panel **E**. CLAUDE.md rule 2: never infer
a panel letter from a directory name. And `03_C/create_tcga_ceacam_scatter.py`
is a *different file with the same name* drawing printed panel B; the two are
kept in their own directories and must not be merged.

Version A is
`03_Revised_Panels/Main_Figures/03_Figure_3/03_D/create_tcga_ceacam_scatter.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every value read, every filter, every statistic and every string is Version A's.
The drawing code is the same code: the 18-gene Tex signature, the weighted-KS
ssGSEA, the z-scoring, the clinical merge, the per-subplot zero-expression
filter, the Spearman test and the P-value formatting are untouched.

  printed panel  Figure 3 E     (PROVENANCE.csv; NOT inferred from "03_D")
  printed rect   47.7 x 21.8 mm   (panel_rects.csv)
  Version B box  107.0 x 57.0 mm

MARK
    Version A drew a 32.0 x 16.0 cm canvas (SCALE = 4) and passed
    `fontscale = SCALE` into `plot_scatter`. Its smallest body type is the
    in-axes legend at `fontsize=4.5 * fontscale` - smaller than the
    5 * fontscale tick labels - so SMALL_PT = 4.5 and, by PANEL_SPEC.md,

        MARK = tick_pt / (SMALL_PT * SCALE) = 7 / 18 = 0.389
        AREA = MARK ** 2                    = 0.151

    Marker area 20*fontscale, marker edge 0.5*fontscale, the dashed regression
    line 0.8*fontscale, the spines 0.5*fontscale, the tick width/length
    0.5*fontscale / 3*fontscale and the two legend keys (markersize 4, edge
    0.5) are each multiplied by it, so all keep their Version A size *relative
    to the type*.

THE LEGEND KEYS ARE DIVIDED BY legend.markerscale
    cnsplots sets `legend.markerscale = 0.5`; Version A ran at matplotlib's
    default 1.0. Carrying the key sizes over unchanged would therefore halve
    them again on top of MARK and close the open 'Alive' ring into a dot. Each
    key is sized so that *after* the scale it prints at MARK times its Version
    A size. Same handles, same labels, same colours. (Same treatment as the
    stage-2 exemplar S8_F.)
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
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import TCGA_BAYESPRISM_EPI, TCGA_BULK_TUMOR_ONLY, TCGA_RAW_DIR, TCGA_CLINICAL
import panel_style_cns as style  # noqa: E402

# ==============================================================================
# Configuration
# ==============================================================================
# Version A's canvas multiplier. Version B draws 1:1, so SCALE survives only to
# reproduce the exact numbers Version A set for its non-type lengths.
SCALE = 4
SMALL_PT = 4.5                       # Version A's smallest body type (legend)

# See MARK above. `style.tick_pt()` reads cnsplots' own setting, so the factor
# is derived, never a literal.
MARK = style.tick_pt() / (SMALL_PT * SCALE)   # 0.389, length multiplier
AREA = MARK ** 2                              # 0.151, area multiplier

PRINTED_MM = (47.7, 21.8)            # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 107.0, 57.0
# Millimetres of paper. left holds the y label + its tick labels, top the
# two-line title, bottom the x label + its tick labels; wspace 1/3 of a column
# leaves 13 mm between the axes for the right subplot's own y furniture.
MARGIN = dict(left=14.0, right=2.0, top=9.0, bottom=9.0, wspace=1.0 / 3.0)

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
                   s=20*fontscale*AREA, alpha=0.85, linewidths=0.5*fontscale*MARK, zorder=3)

    # Regression line (above dots)
    valid = np.isfinite(x) & np.isfinite(y)
    xv, yv = np.array(x)[valid], np.array(y)[valid]
    if len(xv) >= 3 and np.std(xv) > 0:
        slope, intercept = np.polyfit(xv, yv, 1)
        x_line = np.linspace(xv.min(), xv.max(), 100)
        ax.plot(x_line, slope * x_line + intercept, 'k--', linewidth=0.8*fontscale*MARK, alpha=0.6, zorder=5)

    # Stats — 1 sig digit (floor), scientific for very small P
    import math
    _e = math.floor(math.log10(p_val)); _c = int(p_val / 10**_e)
    if _e >= -3:
        p_str = f'P = {_c * 10**_e:.{-_e}f}'
    else:
        p_str = f'P = {_c}' + r'$\times 10^{' + str(_e) + r'}$'

    ax.set_title(f'{title}\nρ = {r_val:.2f}, {p_str}', linespacing=1.4)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for spine in ['bottom', 'left']:
        ax.spines[spine].set_linewidth(0.5*fontscale*MARK)
    ax.tick_params(axis='both', width=0.5*fontscale*MARK, length=3*fontscale*MARK)
    ax.set_box_aspect(1)

    # Compact legend (Dead filled red, Alive open gray). `key` undoes cnsplots'
    # legend.markerscale so each key prints at MARK x its Version A size.
    key = 1.0 / plt.rcParams['legend.markerscale']
    h_dead = mlines.Line2D([], [], marker='o', color='none', markerfacecolor=COLOR_DEAD,
                           markeredgecolor=COLOR_DEAD,
                           markersize=4*fontscale/SCALE*MARK*key,
                           label='Dead')
    h_alive = mlines.Line2D([], [], marker='o', color='none', markerfacecolor='none',
                            markeredgecolor=COLOR_ALIVE,
                            markersize=4*fontscale/SCALE*MARK*key,
                            markeredgewidth=0.5*MARK, label='Alive')
    ax.legend(handles=[h_dead, h_alive], loc='upper left',
              frameon=False, handletextpad=0.3, borderpad=0.2)

    return r_val, p_val


def main():
    family = style.apply()
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

    df = load_data()
    n = len(df)

    # Panel D: CEACAM vs Tex ssGSEA only (filter zero-CEACAM samples per subplot)
    panel = {
        'name': 'tcga_ceacam_tex',
        'combos': [
            ('CEACAM5', 'tex_ssgsea', 'Epi. CEACAM5 (log$_2$+1)', 'Tex ssGSEA (z-score)'),
            ('CEACAM6', 'tex_ssgsea', 'Epi. CEACAM6 (log$_2$+1)', 'Tex ssGSEA (z-score)'),
        ],
    }

    results = []
    fig, axes = style.subplots_mm(PANEL_W_MM, PANEL_H_MM, 1, 2)
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

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, OUTPUT_DIR / panel['name'])
    print(f"Saved: {OUTPUT_DIR / panel['name']}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")

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
