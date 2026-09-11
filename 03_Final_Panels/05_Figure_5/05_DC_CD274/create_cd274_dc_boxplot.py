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
    # The two response groups, named without their sample counts: at
    # 6 pt "(n=5)" sets 5.0 mm and the two ticks are about 3 mm apart,
    # so the counts cannot be printed here without overprinting each
    # other. They are stated in the figure legend instead.
    ax.set_xticklabels(['R', 'NR'])
    print(f"  tick labels R, NR; counts R n={len(r_data)}, "
          f"NR n={len(nr_data)} -> figure legend")
    # Named once, beside the leftmost box of the group.
    ax.set_ylabel('')
    # The significance annotation is placed above the topmost datum,
    # which autoscale does not see; the title is lifted clear of it.
    ax.set_title('DC', pad=style.tick_pt())

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_ylim(ax.get_ylim()[0], bh + 0.15 * y_range)

    style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_CELL,
                      reserve_letter=False)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL) if False else []
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")
    stem = 'cd274_dc_boxplot'
    style.save_panel(fig, OUT_DIR / stem)
    print(f"  Saved: {stem}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
