#!/usr/bin/env python3
"""
Figure 3 panel L, RESTYLED (Version B) - distance to immune across the spatial spots
of one Visium section.

NOTE THE PANEL LETTER. This directory is `03_L` and PROVENANCE.csv
confirms it holds printed panel **L** - but the letter was still looked up
there rather than read off the directory name (CLAUDE.md rule 2).

Version A is
`03_Final_Panels/03_Figure_3/03_L/create_spatial_dist_immune.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every spot, every coordinate, every colour value, every colour map and every
string is Version A's. The drawing code is the same code: the same section, the
same `distance_to_immune` column, the same `cmap='plasma'`, the same `alpha`, the same
`fraction`/`pad` colour bar and the same equal aspect.

  printed panel  Figure 3 L     (PROVENANCE.csv; NOT inferred from "03_L")
  printed rect   33.2 x 26.4 mm   (panel_rects.csv)
  Version B box  62.0 x 48.0 mm

MARK
    Version A drew a 20.0 x 20.0 cm canvas (SCALE = 4). Its smallest body type
    is the axis tick labels and the colour bar tick labels, both at
    `labelsize=6 * SCALE`, so SMALL_PT = 6 and, by PANEL_SPEC.md,

        MARK = tick_pt / (SMALL_PT * SCALE) = 7 / 24 = 0.292
        AREA = MARK ** 2                    = 0.085

    The spot area (`s=6`), the spines, the colour bar outline and both sets of
    tick width/length were literal points on the 4x canvas and are multiplied
    by it, so each keeps its Version A size *relative to the type*.

THE MARGINS ARE SET BEFORE THE COLOUR BAR IS MADE
    `plt.colorbar(ax=ax)` shrinks `ax` by `fraction + pad` and puts the bar in
    the strip it frees. `subplots_adjust` (which is what `style.margins_mm`
    calls) re-derives a subplot axes' position from its gridspec cell, so
    calling it afterwards would restore `ax` to the full cell and the bar would
    sit on top of the data. Version A avoided this by calling `tight_layout()`,
    which Version B cannot use because it re-fits the canvas. The only change
    is the order of two layout calls; nothing plotted moves.
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import SPATIAL_SPOT_DATA
import panel_style_cns as style  # noqa: E402

# Version A's canvas multiplier. Version B draws 1:1, so SCALE survives only to
# reproduce the exact numbers Version A set for its non-type lengths.
SCALE = 4
SMALL_PT = 6.0                       # Version A's smallest body type

# See MARK above. `style.tick_pt()` reads cnsplots' own setting, so the factor
# is derived, never a literal.
MARK = style.tick_pt() / (SMALL_PT * SCALE)   # 0.292, length multiplier
AREA = MARK ** 2                              # 0.085, area multiplier

PRINTED_MM = (33.2, 26.4)          # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 62.0, 48.0
# Millimetres of paper. left holds the y label and its five-digit tick labels,
# right the colour bar, its tick labels and its own label, bottom the x label
# and its tick labels, top the title.
MARGIN = dict(left=13.0, right=13.0, top=5.0, bottom=10.0)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SAMPLE_NAME = 'sample_03'


def main():
    family = style.apply()
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

    print("=" * 60)
    print("Panel I: Spatial Distance to Immune")
    print("=" * 60)

    print("\n[1/2] Loading spot data...")
    df = pd.read_csv(SPATIAL_SPOT_DATA)
    sample_df = df[df['sample'] == SAMPLE_NAME].copy()
    print(f"  {SAMPLE_NAME}: {len(sample_df)} spots")

    print("\n[2/2] Creating visualization...")
    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)
    style.margins_mm(fig, **MARGIN)

    sc = ax.scatter(sample_df['x'], sample_df['y'],
                    c=sample_df['distance_to_immune'],
                    cmap='plasma', s=6 * AREA, alpha=0.8)
    cbar = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.solids.set_rasterized(False)
    cbar.set_label('Dist to Immune')
    cbar.ax.tick_params(width=1.0 * MARK, length=4 * MARK)
    cbar.outline.set_linewidth(1.0 * MARK)

    ax.set_title('Distance to Immune')
    ax.set_xlabel('X coordinate')
    ax.set_ylabel('Y coordinate')
    ax.tick_params(axis='both', width=1.0 * MARK, length=4 * MARK)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0 * MARK)
    ax.set_aspect('equal')

    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    out = os.path.join(BASE_DIR, 'spatial_dist_immune')
    style.save_panel(fig, out)
    print(f"\n  Saved: {out}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")

    print("\nDone!")


if __name__ == '__main__':
    main()
