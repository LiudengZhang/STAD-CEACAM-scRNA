#!/usr/bin/env python3
"""
Figure 4 panel C, RESTYLED (Version B) - IM-MoMac module proportion in
post-treatment stomach samples, responders versus non-responders.

Version A is
`03_Revised_Panels/Main_Figures/04_Figure_4/04_C/create_module2_boxplot.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every value read, every filter, every statistic and every string is Version A's.
The drawing code is the same code.

  printed panel  Figure 4 C     (PROVENANCE.csv; NOT inferred from "04_C")
  printed rect   31.0 x 39.2 mm    (panel_rects.csv)
  Version B box  40.0 x 48.0 mm

MARK
    Version A drew 3.0 x 3.2 cm at SCALE = 4 (12.0 x 12.8 cm of canvas) and
    set every font from `FONT_SIZE = 5 * SCALE`. The assembler fitted the saved
    311.3 x 334.0 pt SVG into 31.0 x 39.2 mm, a fit of 0.2823, so that 20 pt
    type printed at 5.65 pt. So SMALL_PT = 5, SCALE = 4 and

        MARK = tick_pt / (SMALL_PT * SCALE) = 7 / 20 = 0.35
        AREA = MARK ** 2

    Every length Version A set in points is written here exactly as Version A
    wrote it with `* MARK` appended, so the arithmetic is auditable line by
    line: the box, whisker, cap and median widths, the flier marker and its
    edge, and the significance bracket.

THE FRAME IS cnsplots'
    Version A also set `axes.linewidth`, both spine widths and the tick width
    and length from `SCALE`. Those are the axes frame, which cnsplots declares
    (`axes.linewidth` 0.5, tick size 2, width 0.6, pad 1) and `describe()`
    reports; overriding them per panel is the divergence the restyle exists to
    remove. They are dropped and cnsplots' values stand - the same choice
    `_restyled/S7_Cohort_Statistics/S7_B/create_S7_B_forest.py` made in stage 2.
    MARK governs the marks that draw the data, not the ruler they are drawn on.

Input : FIG4_MODULE_PROPORTIONS; FIG4_CLINICAL_METADATA (00_Config/paths.py)
Output: this directory / module2_post_response_boxplot_k5.{svg,pdf,png}
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import mannwhitneyu
import sys

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import *                       # noqa: E402,F401,F403
import panel_style_cns as style           # noqa: E402

# Version A's canvas convention, kept only so MARK can be derived from it.
SCALE = 4
SMALL_PT = 5.0                            # Version A's `FONT_SIZE = 5 * SCALE`

PRINTED_MM = (31.0, 39.2)                 # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 40.0, 48.0
MARGIN = dict(left=10.0, right=2.0, top=2.0, bottom=5.5)

# Standard colors from CLAUDE.md
COLOR_R = '#0072B2'    # Blue for Responder
COLOR_NR = '#D55E00'   # Vermillion for Non-responder
PALETTE = [COLOR_R, COLOR_NR]  # R first (left), NR second (right)

# Darker shades for median lines
MEDIAN_COLOR_R = '#005689'   # Darker blue
MEDIAN_COLOR_NR = '#A34700'  # Darker vermillion

def main():
    print("=" * 60)
    print("Panel C: IM-MoMac Response Boxplot")
    print("=" * 60)

    script_dir = Path(__file__).parent

    family = style.apply()
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.4f}")

    # Load data from central config paths
    print("\nLoading data...")
    module_props = pd.read_csv(FIG4_MODULE_PROPORTIONS, index_col=0)
    clinical_data = pd.read_csv(FIG4_CLINICAL_METADATA)

    if 'Sample ID' in clinical_data.columns:
        clinical_data = clinical_data.set_index('Sample ID')

    merged = module_props.join(clinical_data, how='inner')

    # Filter for post-treatment stomach samples
    filtered = merged[
        (merged['Sample site'] == 'Stomach') &
        (merged['Treatment phase'] == 'Post') &
        (merged['stomach_post_grouping'].isin(['No-response', 'Responsed']))
    ].copy()

    print(f"  Responsed (R): {len(filtered[filtered['stomach_post_grouping'] == 'Responsed'])}")
    print(f"  No-response (NR): {len(filtered[filtered['stomach_post_grouping'] == 'No-response'])}")

    # Statistical test (Mann-Whitney U)
    group_r = filtered[filtered['stomach_post_grouping'] == 'Responsed']['Module_2'].values
    group_nr = filtered[filtered['stomach_post_grouping'] == 'No-response']['Module_2'].values

    stat, p_value = mannwhitneyu(group_r, group_nr, alternative='two-sided')
    print(f"  Mann-Whitney p-value: {p_value:.4f}")

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

    # Create boxplot with R on left, NR on right
    sns.boxplot(
        data=filtered,
        x='stomach_post_grouping',
        y='Module_2',
        order=['Responsed', 'No-response'],  # R first (left)
        palette=PALETTE,
        ax=ax,
        boxprops=dict(edgecolor='black', linewidth=1.0 * SCALE * MARK),
        medianprops=dict(color='black', linewidth=1.5 * SCALE * MARK),  # Will be overridden
        whiskerprops=dict(linewidth=1.0 * SCALE * MARK),
        capprops=dict(linewidth=1.0 * SCALE * MARK),
        flierprops=dict(marker='o', markerfacecolor='white', markeredgecolor='black',
                        markersize=4 * SCALE * MARK, markeredgewidth=1 * SCALE * MARK),
        width=0.6
    )

    # Color median lines to match box colors (darker shades)
    for i, artist in enumerate(ax.patches):
        if i == 0:
            ax.lines[4].set_color(MEDIAN_COLOR_R)  # First median
        elif i == 1:
            ax.lines[9].set_color(MEDIAN_COLOR_NR)  # Second median

    # Update x-tick labels
    ax.set_xticklabels(['R', 'NR'])

    # Add significance bracket
    y_max = filtered['Module_2'].max()
    y_range = filtered['Module_2'].max() - filtered['Module_2'].min()
    bracket_height = y_max + y_range * 0.1
    bracket_top = bracket_height + y_range * 0.05

    ax.plot([0, 0, 1, 1], [bracket_height, bracket_top, bracket_top, bracket_height],
            lw=0.8 * SCALE * MARK, c='black')
    p_text = '***' if p_value < 0.001 else '**' if p_value < 0.01 else '*' if p_value < 0.05 else 'ns'
    ax.text(0.5, bracket_top + y_range * 0.02, p_text,
            ha='center', va='bottom', fontsize=style.tick_pt() * 0.85)

    # Styling
    ax.set_ylabel('Module Proportion')
    ax.set_xlabel('')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_ylim(ax.get_ylim()[0], bracket_top + y_range * 0.15)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, script_dir / 'module2_post_response_boxplot_k5')
    print(f"\nSaved: {script_dir / 'module2_post_response_boxplot_k5'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")

if __name__ == "__main__":
    main()
