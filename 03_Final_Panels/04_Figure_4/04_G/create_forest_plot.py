#!/usr/bin/env python3
"""
Figure 4 panel G - fold change of each MoMac state between responders and
non-responders, with its confidence interval.

  printed panel  Figure 4 G       (PROVENANCE.csv - the directory is "04_G";
                                   do NOT read the directory as the letter)

The source table, the row order, the fold changes, the confidence intervals,
the reference line and the colour rule are unchanged. The canvas and the type
change: the panel is drawn at the millimetre rectangle it prints in and set in
the figure's one type system.

THE ROW NAMES AND THE AXIS TITLE ARE RESTORED
    The printed page names eight distinct rows and titles the axis
    "Fold Change (Post-R/Post-NR)". The label function this script carried cut
    every state name back to its lineage, which printed five rows reading "Mac"
    and two reading "Mono", and the axis title had lost the treatment phase.
    Both are restored to what the page prints; see DISPLAY_NAME. No fold
    change, interval, row position or colour moves - only the strings.

MARK
    The earlier drawing used a canvas four times the printed size and set its
    tick labels from `TICK_FONTSIZE = 5 * SCALE`, so SMALL_PT = 5 and

        MARK = style.tick_pt() / (SMALL_PT * SCALE)
        AREA = MARK ** 2

    Every length the earlier drawing set in points is written here exactly as
    it was written there with `* MARK` appended - the interval line, the caps,
    the reference line, the marker edge and the axis-label pad - and the square
    marker, which is an area, with `* AREA`. Tick widths and lengths and spine
    widths are not scaled: those are the axes frame, which cnsplots declares.

THE DECADES ARE WRITTEN OUT
    The x axis is logarithmic and labels its decades in scientific notation.
    Mathtext draws a superscript at 70% of its base, so the exponent of a
    6 pt tick label prints at 4.2 pt, below the floor this figure is set to.
    The decades are written out instead. The tick positions, the scale, the
    limits and the values do not move; only the notation changes, and every
    substitution is declared in 00_Config/shared/labels.py.

Input : FIG4_MOMAC_FOLDCHANGE (00_Config/paths.py)
Output: this directory / momac_forest_plot_final.{svg,pdf,png}
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from pathlib import Path
import sys

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *                       # noqa: E402,F401,F403
import panel_style_cns as style           # noqa: E402
import slots                              # noqa: E402

SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 5.0                  # the earlier tick type, before * SCALE

PANEL_LETTER = "G"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(4, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(4, PANEL_LETTER)

# Colors (matching Figure 2 Panel A diverging scheme)
COLOR_PROTECTIVE = '#0072B2'  # Blue for OR < 1
COLOR_RISK = '#D55E00'        # Red/Orange for OR > 1
COLOR_TOTAL = '#009E73'       # Green for Total MoMac
COLOR_NEUTRAL = '#666666'     # Gray for non-significant


#: The name the figure prints for each row. The eight states are named here
#: as panels A, D and E of this figure name them, so that one figure names a
#: cell state one way. Cutting a state name back to its lineage instead leaves
#: five rows reading "Mac" and two reading "Mono"; the page prints eight
#: distinct names, and the page is the ground truth.
DISPLAY_NAME = {
    'C0_Mac_Classic_TREM2':          'Mac_TREM2',
    'C1_Mono_Classic_CD14':          'Mono_CD14',
    'C2_MoMac_Intermediate_HLA-DRA': 'MoMac_Inter',
    'C3_Mac_Inflam_IL1B':            'MoMac_IL1B',
    'C4_Mono_Alternative_CD16':      'Mono_CD16',
    'C5_Mac_Prolif_MKI67':           'Mac_Prolif',
    'C6_Mac_Metallothionein_MT1G':   'Mac_MT1G',
    'Total_MoMac':                   'Total MoMac',
}


def get_short_label(cell_state):
    """The name this figure prints for one cell state."""
    try:
        return DISPLAY_NAME[cell_state]
    except KeyError:
        raise KeyError(
            f"{cell_state!r} has no printed name. Every row of this panel is "
            f"named on the published page; a state that is not in DISPLAY_NAME "
            f"is a state the panel has never printed.") from None


def decade(value, _pos):
    """A power of ten written out, so no glyph on the axis is a superscript."""
    return f"{value:g}"


def main():
    print("=" * 60)
    print("Panel G: Forest Plot (Diverging Colors)")
    print("=" * 60)

    script_dir = Path(__file__).parent

    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.4f}")

    # Load data from central config paths
    results_df = pd.read_csv(FIG4_MOMAC_FOLDCHANGE)
    results_df = results_df.sort_values('Fold_Change', ascending=False)

    print(f"\nLoaded {len(results_df)} cell states")

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

    y_positions = np.arange(len(results_df))

    for idx, row in enumerate(results_df.itertuples()):
        is_total = row.Cell_State == 'Total_MoMac'
        fc = row.Fold_Change
        ci_l = row.CI_Lower
        ci_u = row.CI_Upper
        p_val = row.P_Value

        if np.isnan(fc) or np.isnan(ci_l) or np.isnan(ci_u):
            continue

        # Determine color based on fold change direction and Total MoMac
        if is_total:
            color = COLOR_TOTAL  # Green for Total MoMac
        else:
            color = COLOR_PROTECTIVE if fc < 1 else COLOR_RISK

        # Marker and line properties (same for all)
        line_width = 1.0 * SCALE * MARK
        marker_size = 35 * SCALE * AREA

        # Plot whisker (CI)
        ax.plot([ci_l, ci_u], [idx, idx], color=color, linewidth=line_width,
                alpha=0.85, zorder=2, solid_capstyle='butt')

        # Add caps
        cap_height = 0.15
        ax.plot([ci_l, ci_l], [idx - cap_height, idx + cap_height], color=color,
                linewidth=line_width, alpha=0.85, zorder=2)
        ax.plot([ci_u, ci_u], [idx - cap_height, idx + cap_height], color=color,
                linewidth=line_width, alpha=0.85, zorder=2)

        # Plot square marker (same style for all, edge = fill color)
        ax.scatter(fc, idx, s=marker_size, marker='s', color=color, alpha=0.85,
                  edgecolors=color, linewidths=1.2 * SCALE * MARK, zorder=4)

    # Reference line at FC=1.0
    ax.axvline(x=1, color='black', linestyle='--',
               linewidth=1.0 * SCALE * MARK, alpha=0.7, zorder=1)

    # Y-axis labels
    short_labels = [get_short_label(state) for state in results_df['Cell_State']]
    ax.set_yticks(y_positions)
    ax.set_yticklabels(short_labels)

    # X-axis label
    ax.set_xlabel('Fold Change (Post-R/Post-NR)',
                  labelpad=5 * SCALE * MARK)

    # X-axis settings
    ax.set_xscale('log')
    ax.set_xlim(0.04, 12)
    ax.xaxis.set_major_formatter(FuncFormatter(decade))

    # Grid - only horizontal lines, vertical reference line at x=1 is kept via axvline
    ax.grid(axis='y', alpha=0.1, linestyle=':', zorder=0, color='gray')

    # Spines
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
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    stem = 'momac_forest_plot_final'
    style.save_panel(fig, script_dir / stem)
    print(f"\nSaved: {stem}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == "__main__":
    main()
