#!/usr/bin/env python3
# REVISED FOR CIR-26-0753-ET, reviewer 1 point R1.3c:
# the test is two-sided and the annotation reports the exact P value.
# Original one-tailed version: Round_5/03_Final_Panels/02_Figure_2/02_M/create_ihc_combined_boxplot.py
"""
Figure 2, printed panel N, RESTYLED (Version B) - combined CEACAM5+6 IHC
staining, R versus NR.

Two-sided Mann-Whitney U test, n=4 R, n=4 NR (all pre-treatment).
Data from color deconvolution (Ruifrok & Johnston 2001) via skimage rgb2hed.
DAB OD threshold = 0.02. Source: 02_Preparation_for_Panels/IHC/quantify_ceacam_ihc.py

Version A is
`03_Revised_Panels/Main_Figures/02_Figure_2/02_M/create_ihc_combined_boxplot.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md
allows: the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every value read, every filter, every statistic and every string is Version A's.
The drawing code is the same code.

  printed panel  Figure 2 N     (PROVENANCE.csv; NOT inferred from "02_M")
  printed rect   31.7 x 29.4 mm    (panel_rects.csv)
  Version B box  42.0 x 42.0 mm

MARK
    Version A drew at SCALE = 4 (3.2 x 3.0 cm x 4 = 128 x 120 mm) and its
    smallest body type is the 6 pt nominal on the y label, the tick labels and
    the bracket P value - `6 * SCALE`. So

        SCALE = 4, SMALL_PT = 6
        MARK  = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 24 = 0.29167
        AREA  = MARK ** 2

    Areas (the strip-point `s = 20`) take AREA, lengths take MARK. Version A's
    spine linewidths and `tick_params(width=..., length=...)` are NOT carried
    over: those are axes furniture, cnsplots has its own settings for them,
    and following the library rather than rescaling the old numbers is the
    standard-methods rule.

    The panel grew from 31.7 x 29.4 mm to 42 x 42 mm: the two-line x tick
    labels and the y label at 7/8 pt need the room.
"""

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import sys

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))

OUTPUT_DIR = Path(__file__).parent

# Colors matching Figure 2 panels J/K
COLOR_R = '#0072B2'
COLOR_NR = '#D55E00'
MEDIAN_R = '#005689'
MEDIAN_NR = '#A34700'

PRINTED_MM = (31.7, 29.4)           # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 42.0, 42.0
MARGIN = dict(left=10.5, right=2.5, top=5.5, bottom=9.0)

SCALE = 4                           # Version A's canvas multiplier
SMALL_PT = 6.0                      # Version A's smallest body type

# Read IHC data from color deconvolution results
from paths import IHC_COLOR_DECONV_CSV  # noqa: E402
import panel_style_cns as style         # noqa: E402
IHC_CSV = IHC_COLOR_DECONV_CSV
df = pd.read_csv(IHC_CSV)

# Pivot: combined CEACAM5 + CEACAM6 per patient
pivot = df.pivot_table(index=['patient', 'group'], columns='marker',
                       values='staining_pct').reset_index()
pivot['combined'] = pivot['CEACAM5'] + pivot['CEACAM6']

R_VALS = pivot[pivot['group'] == 'R']['combined'].values
NR_VALS = pivot[pivot['group'] == 'NR']['combined'].values


def main():
    family = style.apply()
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.5f}")

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

    # Two-sided Mann-Whitney U
    stat, pval = stats.mannwhitneyu(NR_VALS, R_VALS, alternative='two-sided')
    print(f"R  (n={len(R_VALS)}): mean={R_VALS.mean():.2f}%")
    print(f"NR (n={len(NR_VALS)}): mean={NR_VALS.mean():.2f}%")
    print(f"Fold: {NR_VALS.mean() / R_VALS.mean():.2f}x")
    print(f"Mann-Whitney U (two-sided): U={stat:.1f}, P={pval:.4f}")

    # Boxplot
    bp = ax.boxplot(
        [R_VALS, NR_VALS],
        positions=[1, 2],
        widths=0.6,
        patch_artist=True,
        boxprops=dict(linewidth=0.5 * MARK),
        whiskerprops=dict(color='black', linewidth=0.5 * MARK),
        capprops=dict(color='black', linewidth=0.5 * MARK),
        flierprops=dict(markersize=0),  # hide outliers, show strip points instead
    )

    bp['boxes'][0].set_facecolor(COLOR_R)
    bp['boxes'][0].set_edgecolor('black')
    bp['boxes'][1].set_facecolor(COLOR_NR)
    bp['boxes'][1].set_edgecolor('black')
    bp['medians'][0].set_color(MEDIAN_R)
    bp['medians'][0].set_linewidth(0.8 * MARK)
    bp['medians'][1].set_color(MEDIAN_NR)
    bp['medians'][1].set_linewidth(0.8 * MARK)

    # Strip points (essential with small n)
    rng = np.random.default_rng(42)
    jitter_r = rng.uniform(-0.12, 0.12, len(R_VALS))
    jitter_nr = rng.uniform(-0.12, 0.12, len(NR_VALS))

    ax.scatter(1 + jitter_r, R_VALS, c=COLOR_R, s=20 * AREA, zorder=5,
               edgecolors='black', linewidths=0.3 * MARK)
    ax.scatter(2 + jitter_nr, NR_VALS, c=COLOR_NR, s=20 * AREA, zorder=5,
               edgecolors='black', linewidths=0.3 * MARK)

    # Significance bracket
    y_max = max(np.max(R_VALS), np.max(NR_VALS))
    y_bracket = y_max * 1.15

    ax.plot([1, 1, 2, 2],
            [y_bracket, y_bracket * 1.05, y_bracket * 1.05, y_bracket],
            'k-', linewidth=0.5 * MARK)

    p_text = f'P = {pval:.3f}' if pval >= 0.001 else 'P < 0.001'

    ax.text(1.5, y_bracket * 1.08, p_text, ha='center', va='bottom',
            fontsize=style.tick_pt())

    # Labels
    ax.set_title('IHC (CEACAM5+6)')
    ax.set_ylabel('Staining area (%)')
    ax.set_xticks([1, 2])
    ax.set_xticklabels([f'R\n(n={len(R_VALS)})', f'NR\n(n={len(NR_VALS)})'])
    ax.set_ylim(0, y_max * 1.40)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")

    style.save_panel(fig, OUTPUT_DIR / "ihc_combined_boxplot")
    print(f"\nSaved: {OUTPUT_DIR / 'ihc_combined_boxplot'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
