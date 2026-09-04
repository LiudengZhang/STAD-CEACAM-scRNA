#!/usr/bin/env python3
"""
Figure 4 panel H, RESTYLED (Version B) - external validation of the Mac_IL1B
(Mac_3) proportion, tumour versus normal, in two public gastric cohorts.

Version A is
`03_Revised_Panels/Main_Figures/04_Figure_4/04_H/create_external_validation_boxplots.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every value read, every filter, every statistic and every string is Version A's.
The drawing code is the same code. Note in particular that the y axis label is
still 'Mac_3 Proportion': the *printed* panel reads 'Mac_IL1B Proportion'
(PROVENANCE.csv row 44 records the discrepancy) and correcting it here would be
a content change, which a restyle may not make.

  printed panel  Figure 4 H     (PROVENANCE.csv; NOT inferred from "04_H")
                 ONE printed letter, TWO drawings - both files are kept.
  printed rect   28.6 x 51.2 mm    for the pair stacked (panel_rects.csv),
                 so roughly 28.6 x 25.6 mm each
  Version B box  44.0 x 40.0 mm    each

MARK
    Version A drew each boxplot 3.0 x 3.0 cm at SCALE = 4 and set every font
    from `TICK_FONTSIZE = 5 * SCALE`. The assembler fitted the saved
    311.3 x 311.4 pt SVGs into ~28.6 x 25.6 mm, a fit of 0.2331, so that 20 pt
    type printed at 4.66 pt. So SMALL_PT = 5, SCALE = 4 and

        MARK = tick_pt / (SMALL_PT * SCALE) = 7 / 20 = 0.35
        AREA = MARK ** 2

    Every length Version A set in points is written here exactly as Version A
    wrote it with `* MARK` appended.

THE FRAME IS cnsplots'
    Version A's `axes.linewidth`, spine widths and tick width/length are the
    axes frame, which cnsplots declares and `describe()` reports. They are
    dropped so the whole figure carries one frame weight - the same choice
    `_restyled/S7_Cohort_Statistics/S7_B/create_S7_B_forest.py` made in stage 2.
    The axes title also loses `fontweight='normal'`: cnsplots bolds axis titles
    and following the library rather than the panel's own override is the
    standard-methods rule.

Input : FIG4_EXTERNAL_LINGHUA; FIG4_EXTERNAL_KUMAR (00_Config/paths.py)
Output: this directory / mac3_boxplot_linghua.{svg,pdf,png}
        this directory / mac3_boxplot_kumar.{svg,pdf,png}
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
SMALL_PT = 5.0                            # Version A's `TICK_FONTSIZE = 5 * SCALE`

PRINTED_MM = (28.6, 51.2)                 # published rect for the pair
PANEL_W_MM, PANEL_H_MM = 44.0, 40.0       # each of the two drawings
MARGIN = dict(left=11.0, right=2.0, top=5.0, bottom=5.5)

# Colors for Tumor vs Normal comparison
COLOR_TUMOR = '#E377C2'   # Pink for Tumor
COLOR_NORMAL = '#17BECF'  # Teal for Normal

# Darker shades for median lines
MEDIAN_COLOR_TUMOR = '#B85A9A'   # Darker pink
MEDIAN_COLOR_NORMAL = '#0F8A99'  # Darker teal

def create_boxplot(data, title, output_stem, y_label='Mac_3 Proportion'):
    """Create a single boxplot for T vs N comparison."""

    MARK = style.tick_pt() / (SMALL_PT * SCALE)

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

    # Create boxplot - Tumor on left, Normal on right
    sns.boxplot(
        data=data,
        x='Group',
        y='Mac3_Proportion',
        order=['T', 'N'],
        palette=[COLOR_TUMOR, COLOR_NORMAL],
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
            ax.lines[4].set_color(MEDIAN_COLOR_TUMOR)  # First median
        elif i == 1:
            ax.lines[9].set_color(MEDIAN_COLOR_NORMAL)  # Second median

    # Statistical test
    tumor_vals = data[data['Group'] == 'T']['Mac3_Proportion'].values
    normal_vals = data[data['Group'] == 'N']['Mac3_Proportion'].values

    stat, p_value = mannwhitneyu(tumor_vals, normal_vals, alternative='two-sided')

    # Add significance bracket
    y_max = data['Mac3_Proportion'].max()
    y_range = data['Mac3_Proportion'].max() - data['Mac3_Proportion'].min()
    if y_range == 0:
        y_range = 0.1
    bracket_height = y_max + y_range * 0.1
    bracket_top = bracket_height + y_range * 0.08

    ax.plot([0, 0, 1, 1], [bracket_height, bracket_top, bracket_top, bracket_height],
            lw=0.8 * SCALE * MARK, c='black')

    p_text = '***' if p_value < 0.001 else '**' if p_value < 0.01 else '*' if p_value < 0.05 else 'ns'
    ax.text(0.5, bracket_top + y_range * 0.02, p_text,
            ha='center', va='bottom', fontsize=style.tick_pt() * 0.85)

    # Styling
    ax.set_title(title)
    ax.set_ylabel(y_label)
    ax.set_xlabel('')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_ylim(ax.get_ylim()[0], bracket_top + y_range * 0.25)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, output_stem)
    print(f"  Saved: {output_stem.name}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")
    print(f"    T: n={len(tumor_vals)}, mean={tumor_vals.mean():.4f}")
    print(f"    N: n={len(normal_vals)}, mean={normal_vals.mean():.4f}")
    print(f"    p-value: {p_value:.4e}")

def main():
    print("=" * 60)
    print("Panel G: External Validation Boxplots")
    print("=" * 60)

    script_dir = Path(__file__).parent

    family = style.apply()
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; "
          f"MARK {style.tick_pt() / (SMALL_PT * SCALE):.4f}")

    # === Linghua Dataset ===
    print("\n[1/2] Linghua (GSE239676)...")
    df_linghua = pd.read_csv(FIG4_EXTERNAL_LINGHUA)

    # Filter to Primary_tumor vs Adjacent_normal only (exclude Peritoneal_carcinomatosis)
    df_linghua_filtered = df_linghua[df_linghua['Tissue'].isin(['Primary_tumor', 'Adjacent_normal'])].copy()
    df_linghua_filtered['Group'] = df_linghua_filtered['Tissue'].map({
        'Primary_tumor': 'T',
        'Adjacent_normal': 'N'
    })

    print(f"  Original samples: {len(df_linghua)}")
    print(f"  After filtering (T vs N only): {len(df_linghua_filtered)}")

    output_linghua = script_dir / 'mac3_boxplot_linghua'
    create_boxplot(df_linghua_filtered, 'GSE239676', output_linghua)

    # === Kumar Dataset ===
    print("\n[2/2] Kumar (GSE183904)...")
    df_kumar = pd.read_csv(FIG4_EXTERNAL_KUMAR)

    # Already has Normal vs Primary_tumor
    df_kumar['Group'] = df_kumar['Tissue'].map({
        'Primary_tumor': 'T',
        'Normal': 'N'
    })

    output_kumar = script_dir / 'mac3_boxplot_kumar'
    create_boxplot(df_kumar, 'GSE183904', output_kumar)

    print("\n" + "=" * 60)
    print("Panel G complete")
    print("=" * 60)

if __name__ == "__main__":
    main()
