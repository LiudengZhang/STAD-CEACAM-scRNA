#!/usr/bin/env python3
"""
Figure 3 panel N, RESTYLED (Version B) - CEACAM vs immune infiltration,
region-level bar chart.

NOTE THE PANEL LETTER. This directory is `03_N` and PROVENANCE.csv confirms it
holds printed panel **N** - but the letter was still looked up there rather
than read off the directory name (CLAUDE.md rule 2).

Version A is
`03_Revised_Panels/Main_Figures/03_Figure_3/03_N/create_immune_recruitment_barplot.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every value read, every ordering, every colour rule, every significance
threshold and every string is Version A's. The drawing code is the same code:
the same sort on `difference`, the same short-name map, the same
red-above-zero / blue-below-zero rule, the same y limits and the same
`***`/`**`/`*` cut-offs at 0.001 / 0.01 / 0.05.

  printed panel  Figure 3 N     (PROVENANCE.csv; NOT inferred from "03_N")
  printed rect   55.6 x 36.0 mm   (panel_rects.csv)
  Version B box  72.0 x 50.0 mm

MARK
    Version A drew a 28.0 x 22.0 cm canvas (SCALE = 4). Its smallest body type
    is the category tick labels at `fontsize=7 * SCALE` (and the matching
    `tick_params(labelsize=7 * SCALE)`), so SMALL_PT = 7 and, by PANEL_SPEC.md,

        MARK = tick_pt / (SMALL_PT * SCALE) = 7 / 28 = 0.25
        AREA = MARK ** 2                    = 0.0625

    The bar edges, the zero rule, the spines and the tick width/length were
    literal points on the 4x canvas and are each multiplied by it, so all keep
    their Version A size *relative to the type*.

TYPE
    The significance stars lose `fontweight='bold'`: cnsplots bolds axis titles
    and panel letters and nothing else (PANEL_SPEC.md, Type). The stars, their
    positions and their thresholds are unchanged.

    The box grew from 55.6 mm because the nine category labels are set at 7 pt
    and rotated 45 degrees: `Neutrophils` alone is ~15 mm long, so ~11 mm of
    depth below the axis is type before any bar is drawn.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import SPATIAL_REGION_COMPARISON
import panel_style_cns as style  # noqa: E402

# Version A's canvas multiplier. Version B draws 1:1, so SCALE survives only to
# reproduce the exact numbers Version A set for its non-type lengths.
SCALE = 4
SMALL_PT = 7.0                       # Version A's smallest body type

# See MARK above. `style.tick_pt()` reads cnsplots' own setting, so the factor
# is derived, never a literal.
MARK = style.tick_pt() / (SMALL_PT * SCALE)   # 0.25, length multiplier
AREA = MARK ** 2                              # 0.0625, area multiplier

PRINTED_MM = (55.6, 36.0)            # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 72.0, 50.0
# Millimetres of paper. left holds the two-line y label and its tick labels,
# bottom the 45-degree category labels, top the title.
MARGIN = dict(left=15.0, right=3.0, top=6.0, bottom=14.0)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))


def main():
    family = style.apply()
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

    print("=" * 60)
    print("Panel K: Immune Recruitment by CEACAM Region")
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
    ax.set_ylabel('Difference in Proportion\n(CEACAM-high - CEACAM-low)')
    ax.set_title('Immune Recruitment by CEACAM Region')
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

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")

    print("\n[3/3] Saving figure...")
    out = os.path.join(BASE_DIR, 'ceacam_immune_region_summary')
    style.save_panel(fig, out)
    print(f"  Saved: {out}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")

    print("\nDone!")


if __name__ == '__main__':
    main()
