#!/usr/bin/env python3
"""
Figure 3 panel N - the difference in immune cell proportion between CEACAM-high
and CEACAM-low regions, one bar per cell type.

The panel is drawn at the millimetre rectangle it prints in, read from
03_Final_Panels/panel_rects.csv through 00_Config/slots.py, so the type size
set here is the type size printed. Margins are measured from the rendered ink
rather than typed, and the panel-letter corner is left clear for the assembler.

  printed panel  Figure 3 N       (PROVENANCE.csv; the directory name agrees,
                                   but the letter was looked up, not inferred)

Every value read, every ordering, every colour rule, every significance
threshold and every string is the earlier drawing's: the same sort on the
difference, the same short-name map, the same red-above-zero and
blue-below-zero rule, the same y limits and the same cut-offs at 0.001, 0.01
and 0.05.

  N's recorded rect sits 0.109 mm below the top of its own title on the
  published page; the table is written to 0.1 mm and the clip is one grid step.
  It is recorded in KNOWN_RECT_ISSUES.md and the rect is used as it stands.

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the category tick labels - at 7 * SCALE. MARK carries
    the non-type point sizes across to the 1:1 canvas:

        MARK = style.tick_pt() / (SMALL_PT * SCALE)

    The bar edges, the zero rule and the tick width and length are scaled by
    it. Spine widths are not: those are style, and cnsplots sets them.

THE AXIS LABEL AND THE TITLE ARE RE-WRAPPED
    Both carry every word they carried before; only where the lines break
    changes. The axis label is rotated, so its length is vertical: on one and
    two lines it sets 35.0 mm against a 36.0 mm panel and cannot clear the
    panel-letter corner at the top of that column. On four lines it sets
    17.4 mm and fits inside the plotting box. The title sets 45.9 mm on one
    line and the plotting box is 37 mm wide; a title is centred on that box, so
    a title wider than it cannot be fitted at all. On two lines it sets
    22.9 mm. Both are declared in RENAMES_FIGURE_3.

TYPE
    The significance stars lose their bold weight: cnsplots bolds axis titles
    and panel letters and nothing else. The stars, their positions and their
    thresholds are unchanged.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import SPATIAL_REGION_COMPARISON
import panel_style_cns as style  # noqa: E402
import slots  # noqa: E402

SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 7.0                  # the earlier smallest body type, before * SCALE

MARK = style.tick_pt() / (SMALL_PT * SCALE)   # length multiplier
AREA = MARK ** 2                              # area multiplier

PANEL_LETTER = "N"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(3, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(3, PANEL_LETTER)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

    print("=" * 60)
    print("Figure 3 panel N: immune recruitment by CEACAM region")
    print("=" * 60)

    print("\n[1/3] Loading region-level results...")
    results_df = pd.read_csv(SPATIAL_REGION_COMPARISON)
    results_df = results_df.sort_values('difference', ascending=False)
    print(f"  Loaded results for {len(results_df)} immune cell types")

    print("\n[2/3] Creating visualization...")
    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

    SHORT_NAMES = {
        'Monocytes/Macrophages': 'Mo/Mac',
        'CD4+ T cells': 'CD4+ T',
        'CD8+ T cells': 'CD8+ T',
        'Plasma cells': 'Plasma',
    }

    df_sorted = results_df.sort_values('difference', ascending=False)
    df_sorted['label'] = df_sorted['immune_cell'].map(lambda x: SHORT_NAMES.get(x, x))
    x_pos = np.arange(len(df_sorted))
    colors = ['#d62728' if d > 0 else '#1f77b4' for d in df_sorted['difference']]

    bars = ax.bar(x_pos, df_sorted['difference'], color=colors, edgecolor='black',
                  alpha=0.8, width=0.7, linewidth=1.0 * MARK)

    ax.axhline(y=0, color='black', linestyle='-', linewidth=1.0 * MARK)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(df_sorted['label'], rotation=45, ha='right')
    ax.set_ylabel('Difference in\nProportion\n(CEACAM-high\n- CEACAM-low)')
    ax.set_title('Immune Recruitment\nby CEACAM Region')
    ax.tick_params(axis='both', width=1.0 * MARK, length=4 * MARK)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0 * MARK)

    y_min = df_sorted['difference'].min() - 0.01
    y_max = df_sorted['difference'].max() + 0.02
    ax.set_ylim(y_min, y_max)

    for i, (_, row) in enumerate(df_sorted.iterrows()):
        pval = row['pvalue']
        sig = ''
        if pval < 0.001:
            sig = '***'
        elif pval < 0.01:
            sig = '**'
        elif pval < 0.05:
            sig = '*'

        diff = row['difference']
        if diff >= 0:
            y_pos_text = diff + 0.002
            va = 'bottom'
        else:
            y_pos_text = diff - 0.002
            va = 'top'
        ax.text(i, y_pos_text, sig, ha='center', va=va)

    style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_CELL)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    print("\n[3/3] Saving figure...")
    out = os.path.join(BASE_DIR, 'ceacam_immune_region_summary')
    style.save_panel(fig, out)
    print(f"  Saved: {out}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")

    print("\nDone!")


if __name__ == '__main__':
    main()
