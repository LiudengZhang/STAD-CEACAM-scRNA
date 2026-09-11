#!/usr/bin/env python3
"""
Figure 3 panel G - the neighbourhood epithelial density across the spatial
spots of one Visium section.

The panel is drawn at the millimetre rectangle it prints in, read from
03_Final_Panels/panel_rects.csv through 00_Config/slots.py, so the type size
set here is the type size printed. Margins are measured from the rendered ink
rather than typed, and the panel-letter corner is left clear for the assembler.

  printed panel  Figure 3 G       (PROVENANCE.csv; the directory name agrees,
                                   but the letter was looked up, not inferred)

Every spot, every coordinate, every colour value, every colour map and every
string is the earlier drawing's: the same section, the same
neighborhood_epi_density column, the same colour map, the same alpha, the same
colour bar and the same equal aspect.

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the axis tick labels and the colour bar tick labels -
    at 6 * SCALE. MARK carries the non-type point sizes across to the 1:1
    canvas:

        MARK = style.tick_pt() / (SMALL_PT * SCALE)

    The spot area, the colour bar outline and both sets of tick width and
    length are scaled by it. Spine widths are not: those are style, and
    cnsplots sets them.

THE COLOUR BAR IS MADE BEFORE THE MARGINS ARE FITTED
    `plt.colorbar(ax=ax)` splits the axes' own grid cell between the plot and
    the bar, so both keep a subplot specification and both move together when
    `fit_margins` changes the margins. Made afterwards, the bar would be placed
    against a box the fit had already left.
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import SPATIAL_SPOT_DATA
import panel_style_cns as style  # noqa: E402
import slots  # noqa: E402

SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 6.0                  # the earlier smallest body type, before * SCALE

MARK = style.tick_pt() / (SMALL_PT * SCALE)   # length multiplier
AREA = MARK ** 2                              # area multiplier

PANEL_LETTER = "G"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(3, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(3, PANEL_LETTER)
# The map fills its slot to the left edge and the title is centred over it,
# so the title's first glyph comes to rest against the panel letter with no
# paper between them. The fit is given the keep-out cell plus a millimetre of
# gutter; the cell itself is what is checked afterwards.
LETTER_GUTTER_MM = 1.0
LETTER_FIT = (LETTER_CELL[0] + LETTER_GUTTER_MM, LETTER_CELL[1])

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SAMPLE_NAME = 'sample_03'


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

    print("=" * 60)
    print("Figure 3 panel G: spatial epithelial density")
    print("=" * 60)

    print("\n[1/2] Loading spot data...")
    df = pd.read_csv(SPATIAL_SPOT_DATA)
    sample_df = df[df['sample'] == SAMPLE_NAME].copy()
    print(f"  {SAMPLE_NAME}: {len(sample_df)} spots")

    print("\n[2/2] Creating visualization...")
    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

    sc = ax.scatter(sample_df['x'], sample_df['y'],
                    c=sample_df['neighborhood_epi_density'],
                    cmap='YlOrRd', s=6 * AREA, alpha=0.8)
    cbar = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.solids.set_rasterized(False)
    # The shipped page prints no label on the colour bar; the quantity is
    # named by the panel title above the map. Declared in
    # REMOVALS_FIGURE_3.
    cbar.set_label('')
    cbar.ax.tick_params(width=1.0 * MARK, length=4 * MARK)
    cbar.outline.set_linewidth(1.0 * MARK)

    ax.set_title('Epithelial Density')
    # The shipped page prints neither axis label nor either ruler on this map:
    # the coordinates are Visium pixel positions within one section and carry
    # nothing the reader reads off them. The two labels are declared in
    # REMOVALS_FIGURE_3; the tick labels are the numeric ruler.
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.tick_params(axis='both', width=1.0 * MARK, length=4 * MARK)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0 * MARK)
    ax.set_aspect('equal')

    style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_FIT)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    out = os.path.join(BASE_DIR, 'spatial_epi_density')
    style.save_panel(fig, out)
    print(f"\n  Saved: {out}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")

    print("\nDone!")


if __name__ == '__main__':
    main()
