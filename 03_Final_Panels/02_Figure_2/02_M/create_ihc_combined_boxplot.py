#!/usr/bin/env python3
# REVISED FOR CIR-26-0753-ET, reviewer 1 point R1.3c:
# the test is two-sided and the annotation reports the exact P value.
# The one-tailed version this replaces is the one used in the preprint,
# https://www.biorxiv.org/content/10.64898/2026.03.05.708917
"""
Figure 2 panel N - combined CEACAM5+6 immunohistochemical staining area,
responders versus non-responders.

Two-sided Mann-Whitney U test over the eight pre-treatment specimens. Staining
area comes from colour deconvolution (Ruifrok & Johnston 2001) through
skimage rgb2hed at a DAB optical-density threshold of 0.02, computed by
upstream/IHC/quantify_ceacam_ihc.py.

The panel is drawn at the millimetre rectangle it prints in, read from
03_Final_Panels/panel_rects.csv through 00_Config/slots.py, so the type size
set here is the type size printed. Margins are measured from the rendered ink
rather than typed, and the panel-letter corner is left clear for the assembler.

  printed panel  Figure 2 N       (PROVENANCE.csv; NOT inferred from "02_M")

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the y axis label, the tick labels and the bracket P
    value - at 6 * SCALE. MARK carries the non-type point sizes across to the
    1:1 canvas:

        MARK = style.tick_pt() / (SMALL_PT * SCALE)

    Every non-type length - the box, whisker, cap, median and bracket line
    widths and the marker sizes - takes it. Tick widths and lengths and spine
    widths do not: those are style, and cnsplots sets them.

Every value read, every filter, every statistic and every string is the earlier
drawing's. The drawing code is the same code.
"""

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import sys

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))

OUTPUT_DIR = Path(__file__).parent

# Colours matching the two boxes of printed panel K
COLOR_R = '#0072B2'
COLOR_NR = '#D55E00'
MEDIAN_R = '#005689'
MEDIAN_NR = '#A34700'

SCALE = 4                           # the earlier canvas multiplier
SMALL_PT = 6.0                      # the earlier smallest body type

# Read IHC data from color deconvolution results
from paths import IHC_COLOR_DECONV_CSV  # noqa: E402
import panel_style_cns as style         # noqa: E402
import slots                            # noqa: E402
from cnsfig.boxes import draw_boxes, bracket, ylim_above, assert_no_points  # noqa: E402

PANEL_LETTER = "N"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(2, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(2, PANEL_LETTER)
IHC_CSV = IHC_COLOR_DECONV_CSV
df = pd.read_csv(IHC_CSV)

# Pivot: combined CEACAM5 + CEACAM6 per patient
pivot = df.pivot_table(index=['patient', 'group'], columns='marker',
                       values='staining_pct').reset_index()
pivot['combined'] = pivot['CEACAM5'] + pivot['CEACAM6']

R_VALS = pivot[pivot['group'] == 'R']['combined'].values
NR_VALS = pivot[pivot['group'] == 'NR']['combined'].values


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
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

    # ONE BOX, NO POINTS (2026-09-16, the author's fifth reading: "box
    # plots should look alike throughout; don't show every point"). The
    # eight patients' strip points (default_rng(42)) are no longer drawn -
    # the author's ruling, a declared departure from the published panel -
    # and the fliers the strip had replaced are back as the family's 0.79 mm
    # open circles; the median is black like every other box. The bracket
    # keeps its vertices (line at 1.15 x y_max, arms 5% of that) and prints
    # through p_text_kw, the paper's one P-label owner ("P = 0.06" for the
    # 0.0571 this panel used to print as "P = 0.057"); its ink 0.4 mm above.
    bp = draw_boxes(ax, [R_VALS, NR_VALS], [1, 2], [COLOR_R, COLOR_NR],
                    width=0.6)

    y_max = max(np.max(R_VALS), np.max(NR_VALS))
    _, p_text, _ = bracket(fig, ax, 1, 2, y_max, y_max, pval, kind="pair",
                           lift=0.15, arm=0.05 * 1.15)

    # Labels
    ax.set_title('IHC (CEACAM5+6)')
    ax.set_ylabel('Staining area (%)')
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['Pre-R', 'Pre-NR'])
    ax.set_ylim(0, y_max * 1.40)
    ylim_above(ax, p_text)
    assert_no_points(ax, bp)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_CELL)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(
            f"ink under the panel letter cell: {intruders}")

    style.save_panel(fig, OUTPUT_DIR / "ihc_combined_boxplot")
    print(f"\nSaved: {OUTPUT_DIR / 'ihc_combined_boxplot'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
