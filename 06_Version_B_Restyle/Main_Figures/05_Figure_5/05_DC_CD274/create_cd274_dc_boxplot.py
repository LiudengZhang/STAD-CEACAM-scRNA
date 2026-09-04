#!/usr/bin/env python3
# REVISED FOR CIR-26-0753-ET, reviewer 1 point R1.3c:
# the test is two-sided and the annotation reports the exact P value.
# Original one-tailed version: Round_5/03_Final_Panels/05_Figure_5/05_DC_CD274/create_cd274_dc_boxplot.py
"""
Figure 5 panel J (one of four), RESTYLED (Version B) - PD-L1 (CD274) boxplot,
DC cells, Post-R vs Post-NR.

Version A is
`03_Revised_Panels/Main_Figures/05_Figure_5/05_DC_CD274/create_cd274_dc_boxplot.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK
(areas by AREA), margins are millimetres, and the save is `style.save_panel`.

Every value read, every filter, every threshold, every statistic and every
string is Version A's. This panel seeds its jitter with
`np.random.default_rng(42)`, not `np.random.seed(0)`; that is left exactly as
Version A has it.

  printed panel  Figure 5 J        (PROVENANCE.csv)
                 printed J is four boxplots side by side: 05_I (MoMac),
                 05_J (Epithelial), 05_K (Fibroblast), 05_DC_CD274 (DC).
  printed rect   86.8 x 28.7 mm    (panel_rects.csv, all four together;
                 ~21.7 x 28.7 mm each)
  Version B box  44.0 x 48.0 mm    (as 05_I: this panel also carries a two-line
                 y axis label)

MARK
    Version A drew 3.5 x 5 cm at SCALE = 4 with its smallest body type at
    `5 * SCALE`, so SMALL_PT = 5 and

        MARK = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 20 = 0.35
        AREA = MARK ** 2                                     = 0.1225

Judgement calls: as 05_I - the bold on a significant P is dropped and its size
emphasis taken from the system (body_pt vs tick_pt); the axes title's explicit
`fontweight='normal'` is dropped so cnsplots' bold axis title applies.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import DC_CELLS_H5AD
import panel_style_cns as style

import numpy as np
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

SCALE = 4                       # Version A's canvas multiplier, for MARK only
SMALL_PT = 5.0                  # Version A's smallest body type, before * SCALE
MIN_CELLS = 20

PRINTED_MM = (86.8, 28.7)
PANEL_W_MM, PANEL_H_MM = 44.0, 48.0
MARGIN = dict(left=13.0, right=3.0, top=10.0, bottom=9.0)

COLOR_R = '#2166AC'
COLOR_NR = '#B2182B'

OUT_DIR = Path(__file__).parent


def main():
    family = style.apply()
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

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

    # Mann-Whitney U (two-sided)
    if len(r_data) >= 2 and len(nr_data) >= 2:
        _, pval = mannwhitneyu(nr_data, r_data, alternative='two-sided')
    else:
        pval = np.nan
    print(f"  P-value (two-sided): {pval:.4f}" if not np.isnan(pval) else "  P-value: N/A")

    # Plot
    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

    bp = ax.boxplot([r_data, nr_data], positions=[0, 1], widths=0.5,
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color='black', linewidth=1.5 * SCALE * MARK),
                    whiskerprops=dict(linewidth=1.0 * SCALE * MARK),
                    capprops=dict(linewidth=1.0 * SCALE * MARK),
                    boxprops=dict(linewidth=1.0 * SCALE * MARK))
    bp['boxes'][0].set_facecolor(COLOR_R)
    bp['boxes'][0].set_alpha(0.6)
    bp['boxes'][1].set_facecolor(COLOR_NR)
    bp['boxes'][1].set_alpha(0.6)

    # Jitter points
    rng = np.random.default_rng(42)
    for k, (data, color) in enumerate(zip([r_data, nr_data], [COLOR_R, COLOR_NR])):
        if len(data) > 0:
            jitter = rng.uniform(-0.08, 0.08, len(data))
            ax.scatter([k] * len(data) + jitter, data, c=color, s=20 * SCALE * AREA,
                       edgecolors='white', linewidths=0.3 * SCALE * MARK, alpha=0.85, zorder=3)

    # Significance bracket
    all_vals = np.concatenate([r_data, nr_data])
    y_max = np.max(all_vals)
    y_range = np.max(all_vals) - np.min(all_vals)
    if y_range == 0:
        y_range = 0.1
    bh = y_max + 0.10 * y_range
    ax.plot([0, 0, 1, 1], [bh - 0.02 * y_range, bh, bh, bh - 0.02 * y_range],
            color='black', linewidth=0.8 * SCALE * MARK)

    p_str = f'P = {pval:.3f}' if pval >= 0.001 else 'P < 0.001'
    is_star = pval < 0.05
    ax.text(0.5, bh + 0.02 * y_range, p_str, ha='center',
            fontsize=style.body_pt() if is_star else style.tick_pt())

    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"R\n(n={len(r_data)})", f"NR\n(n={len(nr_data)})"])
    ax.set_ylabel('PD-L1 (CD274)\nExpression')
    ax.set_title('PD-L1 (CD274)\nDC cells')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_ylim(ax.get_ylim()[0], bh + 0.15 * y_range)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    stem = 'cd274_dc_boxplot'
    style.save_panel(fig, OUT_DIR / stem)
    print(f"  Saved: {stem}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
