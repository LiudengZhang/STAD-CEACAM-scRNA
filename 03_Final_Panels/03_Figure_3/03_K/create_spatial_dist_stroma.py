#!/usr/bin/env python3
"""
Figure 3 panel K - the distance to stroma across the spatial spots of one
Visium section.

The panel is drawn at the millimetre rectangle it prints in, read from
03_Final_Panels/panel_rects.csv through 00_Config/slots.py, so the type size
set here is the type size printed. Margins are measured from the rendered ink
rather than typed, and the panel-letter corner is left clear for the assembler.

  printed panel  Figure 3 K       (PROVENANCE.csv; the directory name agrees,
                                   but the letter was looked up, not inferred)

Every spot, every coordinate, every colour value, every colour map and every
string is the earlier drawing's: the same section, the same
distance_to_stroma column, the same colour map, the same alpha, the same
colour bar and the same equal aspect - except for the unit on the colour bar,
which is the subject of the next section.

THE UNIT ON THE COLOUR BAR
    Changed 2026-09-10. `distance_to_stroma` is in the full-resolution image's
    own pixels, the array unit, and the shipped page printed this bar as
    `0 / 2000 / 4000` with no unit anywhere on it. That was survivable while
    panels I and J beside it printed `Distance (a.u.)`. It stopped being
    survivable on 2026-09-10, when I and J were converted and now print
    `Distance (um)`: an array unit is 0.2874 um, so a reader carrying the
    neighbouring panels' unit across to this one reads every number here as
    3.476 times what it is.

    So the values are converted here too, by the one factor
    00_Config/spatial_scale.py derives from the array geometry - imported, not
    repeated, because a second derivation of the same constant is what
    RULES.md rule 5 exists to stop.

    The unit goes in the panel TITLE and not on the bar. REMOVALS_FIGURE_3's
    own reason for the bar carrying no label is that "the quantity each bar
    measures is named by the panel title above the map"; the unit is part of
    naming the quantity, and putting it where the name already is adds no new
    text object.

    It goes on a SECOND LINE, and that is a measurement rather than a
    preference. Three layouts were drawn and their ink read back off the PDFs:

      published            title x 3.13-25.64 mm, y 1.91-4.38; bar at x 25.40,
                           ticks 0 / 2000 / 4000
      "... Stroma (µm)"    title x 0.60-29.17, y 4.12-6.59; ticks 0 / 1000
      "µm" on the bar      title x 0.60-23.10, y 4.02-6.49; ticks 0 / 1000
      "... Stroma\n(µm)"   title x 3.11-25.62, y 0.52-2.99; bar at x 25.40,
                           ticks 0 / 500 / 1000 / 1500

    A one-line title widens to 28.5 mm in a 30.7 mm panel, so its first glyph
    lands inside the panel-letter keep-out cell and `fit_margins` pushes the
    whole panel down by 2.2 mm to clear it - taking 2.2 mm off the colour bar,
    which then carries two labelled ticks instead of three. Labelling the bar
    itself does the same thing from the other side. The two-line title keeps
    the first line's width, and therefore its x position, within 0.02 mm of the
    published page, keeps the bar where it was, and leaves it four ticks.

    This panel therefore no longer reproduces the colour bar printed in
    00_GROUND_TRUTH/figures/Figure 3.pdf, and is not meant to. See
    `SUPERSEDED.md` beside this script - the marker is not `KNOWN_BROKEN.md`
    and carries no `Retest-by:` date, because there is nothing to retest.
    RULES.md rule 1 sets out the two senses of `reproduces_published = no`.

    The map itself does not move in value: the colours come from a `Normalize`
    over the data's own range, so multiplying every value by one positive
    constant leaves the normalised position of every spot - and therefore every
    colour on the map and on the bar - identical. What changes is the bar's
    tick NUMBERS, the added unit, and the layout the fit derives from them.

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
from spatial_scale import cohort_um_per_unit  # noqa: E402
import panel_style_cns as style  # noqa: E402
import slots  # noqa: E402
from cnsfig.layout import pin_frame_mm  # noqa: E402

SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 6.0                  # the earlier smallest body type, before * SCALE

MARK = style.tick_pt() / (SMALL_PT * SCALE)   # length multiplier
AREA = MARK ** 2                              # area multiplier

PANEL_LETTER = "K"
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
    print("Figure 3 panel K: spatial distance to stroma")
    print("=" * 60)

    print("\n[1/2] Loading spot data...")
    df = pd.read_csv(SPATIAL_SPOT_DATA)
    sample_df = df[df['sample'] == SAMPLE_NAME].copy()
    print(f"  {SAMPLE_NAME}: {len(sample_df)} spots")
    # Array units -> micrometres. The factor is taken off the WHOLE
    # cohort table, not this section's slice, for the same reason
    # panels I and J take it off the whole table: it is one cohort
    # constant, and a per-section factor would make this map's unit
    # differ from the boxplots' by up to 2 per cent.
    um_per_unit = cohort_um_per_unit(df)
    print(f"  distances -> micrometres at {um_per_unit:.9f} um "
          f"per array unit")
    sample_df['distance_to_stroma'] = (
        sample_df['distance_to_stroma'] * um_per_unit)

    print("\n[2/2] Creating visualization...")
    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

    # THE MAP IS A FIELD OF SPOTS, NOT A WASH  (2026-09-11)
    #   At 6 * AREA the markers are wide enough to touch their neighbours, and
    #   the section prints as a smear of colour. The published panel shows
    #   every Visium spot separately, with paper between them, which is what
    #   lets a reader see the tissue holes and the spot grid at all. Area is
    #   the square of diameter, so 2.5 * AREA is a marker about two thirds as
    #   wide. Not one spot, coordinate or value moves.
    sc = ax.scatter(sample_df['x'], sample_df['y'],
                    c=sample_df['distance_to_stroma'],
                    cmap='viridis', s=2.5 * AREA, alpha=0.8)
    cbar = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.solids.set_rasterized(False)
    # The shipped page prints no label on the colour bar; the quantity is
    # named by the panel title above the map. Declared in
    # REMOVALS_FIGURE_3.
    cbar.set_label('')
    cbar.ax.tick_params(width=style.RULE_PT, length=4 * MARK)
    cbar.outline.set_linewidth(style.RULE_PT)

    ax.set_title('Distance to Stroma\n(µm)', fontweight='normal')
    # The shipped page prints neither axis label nor either ruler on this map:
    # the coordinates are Visium pixel positions within one section and carry
    # nothing the reader reads off them. The two labels are declared in
    # REMOVALS_FIGURE_3; the tick labels are the numeric ruler.
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.tick_params(axis='both', width=style.RULE_PT, length=4 * MARK)
    # All four spines, which is how the published maps are framed; cnsplots
    # hides top and right by default and the redraw inherited that, leaving
    # these sections with an L-shaped rule instead of a box.
    for spine in ax.spines.values():
        spine.set_linewidth(style.RULE_PT)
        spine.set_visible(True)
    ax.set_aspect('equal', adjustable='datalim')   # frame = axes box; pinned below

    style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_FIT)
    # ONE FRAME LINE FOR K, L, M AND N  (2026-09-14, evening): the frame top
    # at 117.5 mm and its bottom at 143.5 mm on the page, whatever the slot.
    pin_frame_mm(fig, ax, top_mm=7.5, bottom_mm=PANEL_H_MM - 33.5)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    out = os.path.join(BASE_DIR, 'spatial_dist_stroma')
    style.save_panel(fig, out)
    print(f"\n  Saved: {out}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")

    print("\nDone!")


if __name__ == '__main__':
    main()
