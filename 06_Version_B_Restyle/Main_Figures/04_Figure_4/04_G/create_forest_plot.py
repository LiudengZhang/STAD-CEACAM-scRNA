#!/usr/bin/env python3
"""
Figure 4 panel G, RESTYLED (Version B) - forest plot of the MoMac cell-state
fold changes (responder / non-responder) with bootstrap confidence intervals.

Version A is
`03_Final_Panels/04_Figure_4/04_G/create_forest_plot.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every value read, every sort, every fold change, every confidence interval and
every string is Version A's. The drawing code is the same code.

  printed panel  Figure 4 G     (PROVENANCE.csv; NOT inferred from "04_G")
  printed rect   37.8 x 50.3 mm    (panel_rects.csv)
  Version B box  52.0 x 58.0 mm

MARK
    Version A drew 4.5 x 6.4 cm at SCALE = 4 and set its row labels and tick
    labels from `TICK_FONTSIZE = 5 * SCALE`. The assembler fitted the saved
    483.5 x 694.4 pt SVG into 37.8 x 50.3 mm, a fit of 0.2053, so that 20 pt
    type printed at 4.11 pt. So SMALL_PT = 5, SCALE = 4 and

        MARK = tick_pt / (SMALL_PT * SCALE) = 7 / 20 = 0.35
        AREA = MARK ** 2

    Every length Version A set in points is written here exactly as Version A
    wrote it with `* MARK` (or `* AREA` for the marker area) appended: the
    confidence-interval line, its caps, the square marker and its edge, the
    FC = 1 reference rule, and the x axis label pad.

THE BOX
    Eight rows of state labels at 7 pt need about 24 mm of height before the
    axis furniture; the published 50.3 mm carries them at 58 mm
    with the 7 pt x tick labels and the 8 pt axis label below. The width is set
    by 'Total MoMac' / 'MoMac_Inter' at 7 pt (~15 mm) plus the plotting area.
    The published 0.75 aspect is kept (52/58 = 0.90 including the wider label
    gutter the 7 pt labels need).

THE FRAME IS cnsplots'
    Version A's `axes.linewidth`, spine widths and tick width/length are the
    axes frame, which cnsplots declares and `describe()` reports. They are
    dropped so the whole figure carries one frame weight - the same choice
    `_restyled/S7_Cohort_Statistics/S7_B/create_S7_B_forest.py` made in stage 2.

Input : FIG4_MOMAC_FOLDCHANGE (00_Config/paths.py)
Output: this directory / momac_forest_plot_final.{svg,pdf,png}
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path
import sys

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import *                       # noqa: E402,F401,F403
import panel_style_cns as style           # noqa: E402

# Version A's canvas convention, kept only so MARK can be derived from it.
SCALE = 4
SMALL_PT = 5.0                            # Version A's `TICK_FONTSIZE = 5 * SCALE`

PRINTED_MM = (37.8, 50.3)                 # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 52.0, 58.0
MARGIN = dict(left=15.5, right=3.0, top=2.0, bottom=10.0)

# Colors (matching Figure 2 Panel A diverging scheme)
COLOR_PROTECTIVE = '#0072B2'  # Blue for OR < 1
COLOR_RISK = '#D55E00'        # Red/Orange for OR > 1
COLOR_TOTAL = '#009E73'       # Green for Total MoMac
COLOR_NEUTRAL = '#666666'     # Gray for non-significant

def get_short_label(cell_state):
    """Get shortened label for display."""
    if cell_state == 'Total_MoMac':
        return 'Total MoMac'

    parts = cell_state.split('_')
    if len(parts) >= 3:
        short = f"{parts[0]}_{parts[1]}"
    else:
        short = cell_state

    # Strip the Cx_ prefix from the short label
    if '_' in short:
        prefix = short.split('_')[0]
        if len(prefix) >= 2 and prefix[0] == 'C' and prefix[1:].isdigit():
            short = short[len(prefix) + 1:]

    label_mappings = {
        'MoMac': 'Mac',
        'Plasma': 'PB'
    }

    return label_mappings.get(short, short)

def main():
    print("=" * 60)
    print("Panel F: Forest Plot (Diverging Colors)")
    print("=" * 60)

    script_dir = Path(__file__).parent

    family = style.apply()
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
    ax.axvline(x=1, color='black', linestyle='--', linewidth=1.0 * SCALE * MARK,
               alpha=0.7, zorder=1)

    # Y-axis labels
    short_labels = [get_short_label(state) for state in results_df['Cell_State']]
    ax.set_yticks(y_positions)
    ax.set_yticklabels(short_labels)

    # No bold formatting - keep all labels consistent
    # for i, label in enumerate(ax.get_yticklabels()):
    #     if "Total" in label.get_text():
    #         label.set_fontweight('bold')

    # X-axis label
    ax.set_xlabel('Fold Change (R/NR)', labelpad=5 * SCALE * MARK)

    # X-axis settings
    ax.set_xscale('log')
    ax.set_xlim(0.04, 12)

    # Grid - only horizontal lines, vertical reference line at x=1 is kept via axvline
    ax.grid(axis='y', alpha=0.1, linestyle=':', zorder=0, color='gray')

    # Spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, script_dir / 'momac_forest_plot_final')
    print(f"\nSaved: {script_dir / 'momac_forest_plot_final'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")

if __name__ == "__main__":
    main()
