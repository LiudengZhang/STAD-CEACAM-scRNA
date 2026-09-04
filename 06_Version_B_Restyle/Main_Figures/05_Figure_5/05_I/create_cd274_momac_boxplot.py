#!/usr/bin/env python3
# REVISED FOR CIR-26-0753-ET, reviewer 1 point R1.3c:
# the test is two-sided and the annotation reports the exact P value.
# Original one-tailed version: Round_5/03_Final_Panels/05_Figure_5/05_I/create_cd274_momac_boxplot.py
"""
Figure 5 panel J (one of four), RESTYLED (Version B) - CD274 Post-R vs Post-NR
in Monocytes/Macrophages.

Version A is
`03_Revised_Panels/Main_Figures/05_Figure_5/05_I/create_cd274_momac_boxplot.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK
(areas by AREA), margins are millimetres, and the save is `style.save_panel`.

Every value read, every filter, every threshold, every statistic and every
string is Version A's. The drawing code is the same code. `np.random.seed(0)`
stays exactly where Version A put it - the jitter is decoration, but an
unseeded jitter makes the panel unable to reproduce itself.

  printed panel  Figure 5 J        (PROVENANCE.csv; NOT inferred from "05_I")
                 printed J is four boxplots side by side: 05_I (MoMac),
                 05_J (Epithelial), 05_K (Fibroblast), 05_DC_CD274 (DC).
  printed rect   86.8 x 28.7 mm    (panel_rects.csv, all four together;
                 ~21.7 x 28.7 mm each)
  Version B box  44.0 x 48.0 mm

MARK
    Version A drew 3.5 x 5 cm at SCALE = 4 and set its smallest body type at
    `5 * SCALE`, so SMALL_PT = 5 and

        MARK = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 20 = 0.35
        AREA = MARK ** 2                                     = 0.1225

    Every non-type length written on the 4x canvas (box, whisker, cap and
    median line widths, the bracket, the marker edge) is multiplied by MARK and
    the marker area by AREA, so each keeps its size relative to the type.
    Tick widths/lengths and spine widths are NOT rescaled: those are style, and
    cnsplots sets them (axes.linewidth 0.5, ticks size 2 width 0.6).

Judgement calls, stated plainly:
  - `fontweight='bold'` on the significant-P annotation is dropped, per
    PANEL_SPEC ("cnsplots bolds axis titles and panel letters and nothing
    else"). The size distinction Version A drew between a significant and a
    non-significant P is kept, but taken from the system: body_pt (8) versus
    tick_pt (7) instead of 7*SCALE versus 5*SCALE.
  - The axes title's explicit `fontweight='normal'` is dropped so cnsplots'
    bold axis title applies, for the same reason.
"""
import warnings, numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import MOMAC_H5AD
import panel_style_cns as style

warnings.filterwarnings('ignore')

H5AD = MOMAC_H5AD
OUT_DIR = Path(__file__).parent
SCALE = 4                       # Version A's canvas multiplier, for MARK only
SMALL_PT = 5.0                  # Version A's smallest body type, before * SCALE

PRINTED_MM = (86.8, 28.7)       # published rect, all four of panel J together
PANEL_W_MM, PANEL_H_MM = 44.0, 48.0
MARGIN = dict(left=13.0, right=3.0, top=10.0, bottom=9.0)

BOX_COLORS = {'R': '#2166AC', 'NR': '#B2182B'}
MIN_CELLS = 20

def main():
    # The overlaid points are placed with random jitter. Left unseeded it made
    # this panel the only kind in the figure that could not reproduce itself:
    # two consecutive runs of the unchanged script gave three different SVGs
    # (baseline, run 1 and run 2 all differed). The jitter is decoration - no
    # statistic depends on it - but a panel that redraws differently every time
    # cannot be checked, so it is pinned here.
    np.random.seed(0)

    family = style.apply()
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

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

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)
    bp = ax.boxplot([r_data, nr_data], positions=[0, 1], widths=0.5,
                   patch_artist=True, showfliers=False,
                   medianprops=dict(color='black', linewidth=1.5 * SCALE * MARK),
                   whiskerprops=dict(linewidth=1.0 * SCALE * MARK),
                   capprops=dict(linewidth=1.0 * SCALE * MARK),
                   boxprops=dict(linewidth=1.0 * SCALE * MARK))
    bp['boxes'][0].set_facecolor(BOX_COLORS['R']); bp['boxes'][0].set_alpha(0.6)
    bp['boxes'][1].set_facecolor(BOX_COLORS['NR']); bp['boxes'][1].set_alpha(0.6)

    for k, (data, color) in enumerate(zip([r_data, nr_data], [BOX_COLORS['R'], BOX_COLORS['NR']])):
        jitter = np.random.uniform(-0.08, 0.08, len(data))
        ax.scatter([k] * len(data) + jitter, data, c=color, s=20 * SCALE * AREA,
                  edgecolors='white', linewidths=0.3 * SCALE * MARK, alpha=0.85, zorder=3)

    _, pval = stats.mannwhitneyu(nr_data, r_data, alternative='two-sided')
    p_str = f'P = {pval:.3f}' if pval >= 0.001 else 'P < 0.001'
    is_star = pval < 0.05
    y_max = max(np.max(r_data), np.max(nr_data))
    y_range = y_max - min(np.min(r_data), np.min(nr_data))
    bracket_y = y_max + 0.10 * y_range
    ax.plot([0, 0, 1, 1], [bracket_y - 0.02 * y_range, bracket_y,
            bracket_y, bracket_y - 0.02 * y_range], color='black',
            linewidth=0.8 * SCALE * MARK)
    ax.text(0.5, bracket_y + 0.02 * y_range, p_str, ha='center',
            fontsize=style.body_pt() if is_star else style.tick_pt())

    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"R\n(n={len(r_data)})", f"NR\n(n={len(nr_data)})"])
    ax.set_ylabel('PD-L1 (CD274)\nExpression')
    ax.set_title('PD-L1 (CD274)\nMonocytes/Macrophages')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    stem = 'cd274_momac_boxplot'
    style.save_panel(fig, OUT_DIR / stem)
    print(f"R n={len(r_data)}, NR n={len(nr_data)}, P={pval:.4f}")
    print(f"Saved: {stem}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
