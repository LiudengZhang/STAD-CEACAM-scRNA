#!/usr/bin/env python3
# REVISED FOR CIR-26-0753-ET, reviewer 1 point R1.3c:
# the test is two-sided and the annotation reports the exact P value.
# The one-tailed version this replaces is the one used in the preprint,
# https://www.biorxiv.org/content/10.64898/2026.03.05.708917
"""
Figure 2 panel L, CEACAM5 box - BayesPrism epithelial-deconvolved CEACAM5
expression in the PRJEB25780 bulk cohort, responders versus non-responders.

Printed L is one letter over two drawings, 02_N (CEACAM6) and 02_O (CEACAM5,
this file). Each is drawn at its own printed sub-box, read from
03_Final_Panels/slot_subrects.csv through 00_Config/slots.py, so the type
size set here is the type size printed. Margins are measured from the rendered
ink rather than typed. The printed letter L sits inside the first box, not this
one, so this drawing reserves no corner for it.

  printed panel  Figure 2 L, box 2   (PROVENANCE.csv; NOT inferred from "02_O")

MARK
    This is the one Figure 2 panel that was not drawn on a canvas four times
    the printed size: the earlier drawing set its type at print size on a
    27.0 x 24.5 mm canvas, so SCALE = 1, and its smallest body type is the
    bracket P value at 4 pt. MARK carries the non-type point sizes across:

        MARK = style.tick_pt() / (SMALL_PT * SCALE)

    It is above 1 here for the reason it is below 1 elsewhere - this panel's
    type was already near print size, so holding the marks' size relative to
    the type makes them grow rather than shrink. The violin, box, whisker, cap,
    median and bracket line widths and the jittered point areas take it. Tick
    widths and lengths and spine widths do not: those are style, and cnsplots
    sets them.

Every value read, every filter, every statistic and every string is the earlier
drawing's. The drawing code is the same code.
"""

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import TIGER_BAYESPRISM_EPI, TIGER_META       # noqa: E402
import panel_style_cns as style                          # noqa: E402
import slots                                             # noqa: E402
from cnsfig.boxes import box_xlim, frame_fixed          # noqa: E402

OUTPUT_DIR = Path(__file__).parent

# Standard R vs NR colors
COLOR_R = '#2166AC'
COLOR_NR = '#B2182B'

PANEL_LETTER = "L"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(2, PANEL_LETTER, sub=2)

SCALE = 1                           # the earlier canvas multiplier
SMALL_PT = 4.0                      # the earlier smallest body type


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.5f}")

    # Load BayesPrism deconvolved epithelial expression
    print("Loading BayesPrism epithelial expression...")
    epi_expr = pd.read_csv(TIGER_BAYESPRISM_EPI, sep='\t', index_col=0)
    meta_df = pd.read_csv(TIGER_META, sep='\t')
    meta_df = meta_df[meta_df['Treatment'] != 'Normal']

    # Match samples
    common = [s for s in meta_df['sample_id'] if s in epi_expr.index]
    meta_df = meta_df[meta_df['sample_id'].isin(common)].set_index('sample_id')
    epi_expr = epi_expr.loc[common]

    gene = 'CEACAM5'
    r_vals = np.log2(epi_expr.loc[meta_df['response_NR'] == 'R', gene].values + 1)
    nr_vals = np.log2(epi_expr.loc[meta_df['response_NR'] == 'N', gene].values + 1)

    print(f"R (n={len(r_vals)}): mean={np.mean(r_vals):.3f}")
    print(f"NR (n={len(nr_vals)}): mean={np.mean(nr_vals):.3f}")

    stat, pval = stats.mannwhitneyu(nr_vals, r_vals, alternative='two-sided')
    print(f"Mann-Whitney U (two-sided): U={stat:.0f}, P={pval:.4f}")

    # Plot
    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

    data = [r_vals, nr_vals]
    positions = [1, 2]
    colors = [COLOR_R, COLOR_NR]

    # Violin
    vp = ax.violinplot(data, positions=positions, widths=0.7,
                       showmeans=False, showmedians=False, showextrema=False)
    for i, body in enumerate(vp['bodies']):
        body.set_facecolor(colors[i])
        body.set_edgecolor('black')
        body.set_linewidth(style.EDGE_PT)
        body.set_alpha(0.7)

    # Boxplot
    bp = ax.boxplot(data, positions=positions, widths=0.15, patch_artist=True,
                    showfliers=False,
                    boxprops=dict(facecolor='white', linewidth=style.RULE_PT),
                    whiskerprops=dict(color='black', linewidth=style.RULE_PT),
                    capprops=dict(color='black', linewidth=style.RULE_PT),
                    medianprops=dict(color='black', linewidth=style.RULE_PT))

    # Jittered points
    np.random.seed(42)
    for i, (pos, vals, color) in enumerate(zip(positions, data, colors)):
        jitter = np.random.uniform(-0.08, 0.08, len(vals))
        ax.scatter(pos + jitter, vals, c=color,
                   s=(0.76 * style.PT_PER_MM) ** 2,  # published 0.76 mm across
                   alpha=0.8,
                   edgecolors='white', linewidths=style.EDGE_PT, zorder=3)

    # Significance bracket
    y_max = max(np.max(r_vals), np.max(nr_vals))
    y_bracket = y_max * 1.1
    ax.plot([1, 1, 2, 2],
            [y_bracket, y_bracket * 1.03, y_bracket * 1.03, y_bracket],
            'k-', linewidth=style.RULE_PT)
    p_text = style.p_label(pval)
    ax.text(1.5, y_bracket * 1.05, p_text, ha='center', va='bottom',
            fontsize=style.tick_pt())

    # ONE FRAME FOR THE FOUR BOXES OF K AND L  (2026-09-14, evening)
    #   The author asked for K and L to align with each other and with J.
    #   The four boxes name the same margins in millimetres (cnsfig.boxes
    #   frame_fixed: left 9.0 mm for a box with a y label, 5.4 mm without;
    #   3.8 mm below for the tick labels; one title line above), so their
    #   frames print at the same x, the pairs at the same y, and the x
    #   limits come from box_xlim so the boxes stand off the spines.
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['Pre-R', 'Pre-NR'])
    ax.set_xlim(*box_xlim([1, 2], 0.7, clear=0.25))   # 14 mm frame: the two group names need the pitch
    ax.set_ylim(0, y_max * 1.70)   # headroom: the P string clears the title
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # No y label: the left box of this pair carries it (2026-09-14).
    frame_fixed(fig, ax, title='*CEACAM5* (Epi)', ylabel=None,
                left_mm=5.4, right_mm=0.6, bottom_mm=3.8,
                panel_w_mm=PANEL_W_MM, panel_h_mm=PANEL_H_MM)

    style.save_panel(fig, OUTPUT_DIR / "ceacam5_prjeb25780_boxplot")
    print(f"Saved: {OUTPUT_DIR / 'ceacam5_prjeb25780_boxplot'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
