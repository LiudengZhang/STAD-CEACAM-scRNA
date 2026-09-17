#!/usr/bin/env python3
"""
Figure 5 panel J, box 4 of four - CD274 (PD-L1) expression in dendritic cells,
responders versus non-responders after treatment, one point per sample.

The four boxes print side by side under one panel letter. Each is drawn at its
own printed sub-box, read from 03_Final_Panels/slot_subrects.csv through
00_Config/slots.py - not at a quarter of the panel rect, which is wider than
the four boxes together.

  printed panel  Figure 5 J, box 4   (PROVENANCE.csv - the directory is
                                   "05_DC_CD274"; do NOT read the directory as
                                   the letter)

The h5ad, the stomach/post filter, the R-versus-NR mapping, the twenty-cell
minimum, the sample-level mean, the two-sided Mann-Whitney test and the
bracket geometry are unchanged. `np.random.seed(0)` stays exactly where it was.

Each box keeps its own y tick labels: the four have different ranges, and one
shared limit would redraw the data. The axis is named once, beside the leftmost
box. The box title carries the cell type alone, in the short form the rest of
the figure uses; the gene is named in the shared y-axis label. Both changes are
declared in 00_Config/shared/labels.py.

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type at 5 * SCALE, so

        MARK = style.tick_pt() / (SMALL_PT * SCALE)
        AREA = MARK ** 2

    Box, whisker, cap and median line widths, the bracket width and the jitter
    marker area and edge width are scaled by those. Tick widths and lengths and
    spine widths are not: those are style, and cnsplots sets them.

Judgement calls, stated plainly:
  - `fontweight='bold'` on the significant-P annotation is dropped; cnsplots
    bolds axis titles and panel letters and nothing else. The size distinction
    between a significant and a non-significant P is kept, taken from the
    system: body_pt versus tick_pt.
  - The axes title's explicit `fontweight='normal'` is dropped so cnsplots'
    bold axis title applies, for the same reason.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import DC_CELLS_H5AD
import panel_style_cns as style
import slots
from cnsfig.boxes import finish_two_group, draw_boxes, bracket, ylim_above, assert_no_points

import numpy as np
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 5.0                  # the earlier smallest body type, before * SCALE
MIN_CELLS = 20

PANEL_LETTER = "J"
SUB = 4
PANEL_W_MM, PANEL_H_MM = slots.size_mm(5, PANEL_LETTER, sub=SUB)

# The panel letter is measured from the corner of the whole panel rect, which
# starts above and to the left of this box; only the part of the keep-out that
# falls inside this box has to be reserved.
_rect = slots.rect_mm(5, PANEL_LETTER)
_sub = slots.rect_mm(5, PANEL_LETTER, sub=SUB)
_cw, _ch = slots.letter_cell_mm(5, PANEL_LETTER)
LETTER_CELL = (max(0.0, _cw - (_sub[0] - _rect[0])),
               max(0.0, _ch - (_sub[1] - _rect[1])))

COLOR_R = '#2166AC'
COLOR_NR = '#B2182B'

OUT_DIR = Path(__file__).parent


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
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

    # ONE BOX, NO POINTS (2026-09-16, the author's fifth reading: "box
    # plots should look alike throughout; don't show every point"). The
    # box is cnsfig.boxes.draw_boxes (the 2H/I box, 0.79 mm open fliers);
    # the jittered points (default_rng(42)) are no longer drawn - the
    # author's ruling, a declared departure from the published panel. The
    # bracket keeps its vertices (top line at y_max + 0.10 range, arms 0.02
    # down); the star's ink sits 0.4 mm above the line.
    bp = draw_boxes(ax, [r_data, nr_data], [0, 1], [COLOR_R, COLOR_NR],
                    width=0.5)
    bp['boxes'][0].set_alpha(0.6)
    bp['boxes'][1].set_alpha(0.6)

    # Significance bracket
    all_vals = np.concatenate([r_data, nr_data])
    y_max = np.max(all_vals)
    y_range = np.max(all_vals) - np.min(all_vals)
    if y_range == 0:
        y_range = 0.1
    bh = y_max + 0.10 * y_range
    _, p_text, _ = bracket(fig, ax, 0, 1, y_max, y_range, pval, kind="pair",
                           lift=0.08, arm=0.02)
    # Headroom for the P string (2026-09-14, evening): the frame is pinned
    # to the row's line, so the string must be inside the y limits. 0.24 ->
    # 0.42 on 2026-09-15: at 0.24 the star stood against the title's second
    # line; the bracket and star now sit lower in the frame.
    ax.set_ylim(top=bh + 0.42 * y_range)
    ylim_above(ax, p_text)
    assert_no_points(ax, bp)

    # ONE FAMILY, ONE FRAME LINE  (2026-09-14, evening)
    #   The four boxes of J and the box of L are framed by
    #   cnsfig.boxes.finish_two_group: the tick labels with their counts, the
    #   y label as a rich run at the tick size, the title on the canvas in a
    #   two-line band, the frame's top at 139.8 mm and its bottom at 158.7
    #   mm on the page (sweep_pages.ROW_ALIGN), the boxes off the spines.
    finish_two_group(
        fig, ax, title='DC cells',
        tick_labels=[f'R\n(n={len(r_data)})', f'NR\n(n={len(nr_data)})'],
        positions=[0, 1], width=0.5, ylabel_markup='',
        letter_cell=LETTER_CELL, panel_w_mm=PANEL_W_MM, panel_h_mm=PANEL_H_MM,
        top_extra_mm=-0.34, bottom_mm=6.30)
    print(f"  tick labels R, NR; counts R n={len(r_data)}, "
          f"NR n={len(nr_data)} -> figure legend")
    stem = 'cd274_dc_boxplot'
    style.save_panel(fig, OUT_DIR / stem)
    print(f"  Saved: {stem}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
