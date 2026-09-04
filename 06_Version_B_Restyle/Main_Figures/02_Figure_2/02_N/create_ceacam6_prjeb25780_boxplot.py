#!/usr/bin/env python3
# REVISED FOR CIR-26-0753-ET, reviewer 1 point R1.3c:
# the test is two-sided and the annotation reports the exact P value.
# Original one-tailed version: Round_5/03_Final_Panels/02_Figure_2/02_N/create_ceacam6_prjeb25780_boxplot.py
"""
Figure 2, printed panel L (CEACAM6 half), RESTYLED (Version B) - PRJEB25780
BayesPrism epithelial-deconvolved expression, R versus NR.

Version A is
`03_Final_Panels/02_Figure_2/02_N/create_ceacam6_prjeb25780_boxplot.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md
allows: the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every value read, every filter, every statistic and every string is Version A's.
The drawing code is the same code.

  printed panel  Figure 2 L     (PROVENANCE.csv; NOT inferred from "02_N")
                 Printed L is ONE letter over TWO drawings: 02_N (CEACAM6,
                 this file) and 02_O (CEACAM5). They stay two files, as in
                 Version A.
  printed rect   68.3 x 20.4 mm for the pair (panel_rects.csv), so roughly
                 34 x 20.4 mm for this drawing alone
  Version B box  40.0 x 38.0 mm

MARK
    This is the one Figure 2 panel that was NOT drawn at 4x: Version A calls
    `use_panel_style(font_pt=6, scale=1, ...)` on a 27.0 x 24.5 mm canvas, so
    SCALE = 1. Its smallest body type is the bracket P value at `fontsize=4`.

        SCALE = 1, SMALL_PT = 4
        MARK  = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 4 = 1.75
        AREA  = MARK ** 2

    MARK is above 1 here for exactly the reason it is 0.29 elsewhere: this
    panel's type was already near print size, so holding the marks' size
    relative to the type means the marks grow with it rather than shrink.

    Version A's rcParams block (spine and tick widths, tick sizes) is NOT
    carried over: that is axes furniture, cnsplots has its own settings for it,
    and following the library rather than rescaling the old numbers is the
    standard-methods rule. The title's `pad=2` goes the same way - cnsplots
    sets axes.titlepad.
"""

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import TIGER_BAYESPRISM_EPI, TIGER_META       # noqa: E402
import panel_style_cns as style                          # noqa: E402

OUTPUT_DIR = Path(__file__).parent

# Standard R vs NR colors
COLOR_R = '#2166AC'
COLOR_NR = '#B2182B'

PRINTED_MM = (68.3, 20.4)           # published rect of printed L (both halves)
PANEL_W_MM, PANEL_H_MM = 40.0, 38.0
MARGIN = dict(left=11.0, right=2.5, top=5.5, bottom=6.0)

SCALE = 1                           # Version A's canvas multiplier (not 4 here)
SMALL_PT = 4.0                      # Version A's smallest body type


def main():
    family = style.apply()
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

    gene = 'CEACAM6'
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
        body.set_linewidth(0.5 * MARK)
        body.set_alpha(0.7)

    # Boxplot
    bp = ax.boxplot(data, positions=positions, widths=0.15, patch_artist=True,
                    showfliers=False,
                    boxprops=dict(facecolor='white', linewidth=0.5 * MARK),
                    whiskerprops=dict(color='black', linewidth=0.5 * MARK),
                    capprops=dict(color='black', linewidth=0.5 * MARK),
                    medianprops=dict(color='black', linewidth=0.8 * MARK))

    # Jittered points
    np.random.seed(42)
    for i, (pos, vals, color) in enumerate(zip(positions, data, colors)):
        jitter = np.random.uniform(-0.08, 0.08, len(vals))
        ax.scatter(pos + jitter, vals, c=color, s=6 * AREA, alpha=0.8,
                   edgecolors='white', linewidths=0.3 * MARK, zorder=3)

    # Significance bracket
    y_max = max(np.max(r_vals), np.max(nr_vals))
    y_bracket = y_max * 1.1
    ax.plot([1, 1, 2, 2],
            [y_bracket, y_bracket * 1.03, y_bracket * 1.03, y_bracket],
            'k-', linewidth=0.5 * MARK)
    p_text = f'P = {pval:.3f}' if pval >= 0.001 else 'P < 0.001'
    ax.text(1.5, y_bracket * 1.05, p_text, ha='center', va='bottom',
            fontsize=style.tick_pt())

    # Labels
    ax.set_title(r'$\it{CEACAM6}$ (Epi)')
    ax.set_ylabel('Expression (log2)')
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['R', 'NR'])
    ax.set_ylim(0, y_max * 1.35)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")

    style.save_panel(fig, OUTPUT_DIR / "ceacam6_prjeb25780_boxplot")
    print(f"Saved: {OUTPUT_DIR / 'ceacam6_prjeb25780_boxplot'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
