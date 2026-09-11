#!/usr/bin/env python3
"""
Figure 4 panel C - IM-MoMac module proportion in post-treatment stomach
samples, responders against non-responders.

  printed panel  Figure 4 C       (PROVENANCE.csv - the directory is "04_C";
                                   do NOT read the directory as the letter)

The two source tables, the stomach/post filter, the responder mapping, the
two-sided Mann-Whitney test, the box geometry and the bracket are unchanged.
Only the canvas and the type change: the panel is drawn at the millimetre
rectangle it prints in and set in the figure's one type system.

THE GROUP NAMES AND THE PANEL TITLE ARE RESTORED
    The page names the two groups Post-R and Post-NR and titles the panel
    IM-MoMac. The script named the groups R and NR and set no title. Both
    are restored to what the page prints; the two boxes, the samples behind
    them and the test do not move.

MARK
    The earlier drawing used a canvas four times the printed size and set every
    string from `FONT_SIZE = 5 * SCALE`, so SMALL_PT = 5 and

        MARK = style.tick_pt() / (SMALL_PT * SCALE)

    Every length the earlier drawing set in points - the box, whisker, cap and
    median widths, the flier marker and its edge, and the bracket - is written
    here exactly as it was written there with `* MARK` appended, so the
    arithmetic is auditable line by line. Tick widths and lengths and spine
    widths are not scaled: those are the axes frame, which cnsplots declares.

Judgement call, stated plainly:
  - the significance annotation was set at 0.85 of the body size. At the type
    spec that is 5.1 pt, below the floor the figure is set to, so the factor is
    dropped and the annotation takes its size from the system: body_pt when it
    marks a significant difference, tick_pt when it reads "ns".

Input : FIG4_MODULE_PROPORTIONS; FIG4_CLINICAL_METADATA (00_Config/paths.py)
Output: this directory / module2_post_response_boxplot_k5.{svg,pdf,png}
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import mannwhitneyu
import sys

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *                       # noqa: E402,F401,F403
import panel_style_cns as style           # noqa: E402
import slots                              # noqa: E402

SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 5.0                  # the earlier body type, before * SCALE

PANEL_LETTER = "C"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(4, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(4, PANEL_LETTER)

# Standard colors
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

    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
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
    ax.set_xticklabels(['Post-R', 'Post-NR'])

    # Add significance bracket
    y_max = filtered['Module_2'].max()
    y_range = filtered['Module_2'].max() - filtered['Module_2'].min()
    bracket_height = y_max + y_range * 0.1
    bracket_top = bracket_height + y_range * 0.05

    ax.plot([0, 0, 1, 1], [bracket_height, bracket_top, bracket_top, bracket_height],
            lw=0.8 * SCALE * MARK, c='black')
    p_text = '***' if p_value < 0.001 else '**' if p_value < 0.01 else '*' if p_value < 0.05 else 'ns'
    ax.text(0.5, bracket_top + y_range * 0.02, p_text,
            ha='center', va='bottom',
            fontsize=style.body_pt() if p_value < 0.05 else style.tick_pt())

    # Styling
    ax.set_title('IM-MoMac')
    ax.set_ylabel('Module Proportion')
    ax.set_xlabel('')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_ylim(ax.get_ylim()[0], bracket_top + y_range * 0.15)

    style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_CELL)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    stem = 'module2_post_response_boxplot_k5'
    style.save_panel(fig, script_dir / stem)
    print(f"\nSaved: {stem}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == "__main__":
    main()
