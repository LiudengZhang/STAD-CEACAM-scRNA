#!/usr/bin/env python3
"""Figure 4G: fold change (Post-R / Post-NR) of each MoMac state with its
confidence interval, on a log axis.
Log-axis ticks are written as plain numbers (0.1, 1, 10).

Input:  FIG4_MOMAC_FOLDCHANGE (00_Config/paths.py)
Output: momac_forest_plot_final.{svg,pdf,png} in this directory
Run:    python create_forest_plot.py
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.ticker import FuncFormatter
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *                       # noqa: E402,F401,F403
from cnsfig import Panel, style  # noqa: E402
from cnsfig.layout import pin_frame_mm    # noqa: E402

SCALE = 4                       # marker-size scaling constants
SMALL_PT = 5.0

PANEL = Panel(4, "G", __file__, small_pt=SMALL_PT)

# Fold change below 1 blue, above 1 orange; Total MoMac green
COLOR_PROTECTIVE = '#0072B2'  # Blue for OR < 1
COLOR_RISK = '#D55E00'        # Red/Orange for OR > 1
COLOR_TOTAL = '#009E73'       # Green for Total MoMac


#: Row names, as in panels A, D and E.
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
            f"{cell_state!r} has no display name. Every row of this panel is "
            f"named; add the state to DISPLAY_NAME "
            f"before drawing it.") from None


def decade(value, _pos):
    """A power of ten written out, so no glyph on the axis is a superscript."""
    return f"{value:g}"


def main():
    PANEL.apply()
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2

    results_df = pd.read_csv(FIG4_MOMAC_FOLDCHANGE)
    results_df = results_df.sort_values('Fold_Change', ascending=False)

    fig, ax = style.subplots_mm(PANEL.w_mm, PANEL.h_mm)

    y_positions = np.arange(len(results_df))

    for idx, row in enumerate(results_df.itertuples()):
        is_total = row.Cell_State == 'Total_MoMac'
        fc = row.Fold_Change
        ci_l = row.CI_Lower
        ci_u = row.CI_Upper

        if np.isnan(fc) or np.isnan(ci_l) or np.isnan(ci_u):
            continue

        if is_total:
            color = COLOR_TOTAL  # Green for Total MoMac
        else:
            color = COLOR_PROTECTIVE if fc < 1 else COLOR_RISK

        marker_size = 35 * SCALE * AREA

        ax.plot([ci_l, ci_u], [idx, idx], color=color, linewidth=style.RULE_PT,
                alpha=0.85, zorder=2, solid_capstyle='butt')

        cap_height = 0.15
        ax.plot([ci_l, ci_l], [idx - cap_height, idx + cap_height], color=color,
                linewidth=style.RULE_PT, alpha=0.85, zorder=2)
        ax.plot([ci_u, ci_u], [idx - cap_height, idx + cap_height], color=color,
                linewidth=style.RULE_PT, alpha=0.85, zorder=2)

        ax.scatter(fc, idx, s=marker_size, marker='s', color=color, alpha=0.85,
                  edgecolors=color, linewidths=style.EDGE_PT, zorder=4)

    ax.axvline(x=1, color='black', linestyle='--',
               linewidth=style.RULE_PT, alpha=0.7, zorder=1)

    short_labels = [get_short_label(state) for state in results_df['Cell_State']]
    ax.set_yticks(y_positions)
    ax.set_yticklabels(short_labels)

    ax.set_xlabel('Fold Change\n(Post-R/Post-NR)',
                  labelpad=5 * SCALE * MARK)

    ax.set_xscale('log')
    ax.set_xlim(0.04, 12)
    ax.xaxis.set_major_formatter(FuncFormatter(decade))

    ax.grid(axis='y', alpha=0.1, linestyle=':', zorder=0, color='gray')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    PANEL.fit(fig)
    # Align the plot frame with panel H in the same row.
    pin_frame_mm(fig, ax, top_mm=5.0, bottom_mm=10.5)

    PANEL.finish(fig, 'momac_forest_plot_final')


if __name__ == "__main__":
    main()
