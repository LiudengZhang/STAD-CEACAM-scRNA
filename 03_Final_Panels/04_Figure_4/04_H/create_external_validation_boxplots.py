#!/usr/bin/env python3
"""
Figure 4 panel H - Mac_3 proportion in tumour against normal tissue in two
external cohorts, GSE239676 and GSE183904.

  printed panel  Figure 4 H       (PROVENANCE.csv - the directory is "04_H";
                                   do NOT read the directory as the letter)

The panel prints as two boxes stacked under one letter, so each is drawn at its
own printed sub-box, read from 03_Final_Panels/slot_subrects.csv through
00_Config/slots.py - not at half the panel rect, which is taller than the two
boxes together. Box 1 is GSE239676, box 2 is GSE183904, the order
PROVENANCE.csv names and the order the page prints.

The two source tables, the tissue filters, the T/N mapping, the two-sided
Mann-Whitney test, the box geometry and the bracket are unchanged. Only the
canvas and the type change.

THE AXIS TITLE IS RESTORED, AND SET AS TWO LINES
    The page titles both axes MoMac_IL1B Proportion; the script wrote
    Mac_3 Proportion. The proportions plotted are the same column of the same
    two tables. At the figure's body size the restored title is 25.9 mm of
    rotated type against a 23.2 mm box and does not fit on one line - the fit
    refuses it - so it is set as two. Same words, same axis; only the break is
    new, and it is declared with the restoration.

MARK
    The earlier drawing used a canvas four times the printed size and set its
    tick labels from `TICK_FONTSIZE = 5 * SCALE`, so SMALL_PT = 5 and

        MARK = style.tick_pt() / (SMALL_PT * SCALE)

    Every length the earlier drawing set in points is written here exactly as
    it was written there with `* MARK` appended. Tick widths and lengths and
    spine widths are not scaled: those are the axes frame, which cnsplots
    declares.

Judgement calls, stated plainly:
  - the significance annotation was set at 0.85 of the tick size. At the type
    spec that is 5.1 pt, below the floor the figure is set to, so the factor is
    dropped and the annotation takes its size from the system: body_pt when it
    marks a significant difference, tick_pt when it reads "ns".
  - the panel letter is drawn by the assembler at the top-left of the whole
    panel rect, which starts 4.7 mm to the left of box 1. The part of that
    keep-out that falls inside box 1 is therefore empty, and box 2 is far below
    it; box 1 reserves what remains and box 2 reserves nothing.

Input : FIG4_EXTERNAL_LINGHUA; FIG4_EXTERNAL_KUMAR (00_Config/paths.py)
Output: this directory / mac3_boxplot_linghua.{svg,pdf,png}
        this directory / mac3_boxplot_kumar.{svg,pdf,png}
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
SMALL_PT = 5.0                  # the earlier tick type, before * SCALE

PANEL_LETTER = "H"

# The panel letter is measured from the corner of the whole panel rect, which
# starts above and to the left of both boxes; only the part of the keep-out
# that falls inside a box has to be reserved by it.
_rect = slots.rect_mm(4, PANEL_LETTER)
_cw, _ch = slots.letter_cell_mm(4, PANEL_LETTER)


def letter_cell(sub):
    x0, y0, _x1, _y1 = slots.rect_mm(4, PANEL_LETTER, sub=sub)
    return (max(0.0, _cw - (x0 - _rect[0])), max(0.0, _ch - (y0 - _rect[1])))


# Colors for Tumor vs Normal comparison
COLOR_TUMOR = '#E377C2'   # Pink for Tumor
COLOR_NORMAL = '#17BECF'  # Teal for Normal

# Darker shades for median lines
MEDIAN_COLOR_TUMOR = '#B85A9A'   # Darker pink
MEDIAN_COLOR_NORMAL = '#0F8A99'  # Darker teal


def create_boxplot(data, title, stem, sub, reserve_letter,
                   y_label='MoMac_IL1B\nProportion'):
    """Create a single boxplot for T vs N comparison, at its printed sub-box."""

    w_mm, h_mm = slots.size_mm(4, PANEL_LETTER, sub=sub)
    cell = letter_cell(sub)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)

    fig, ax = style.subplots_mm(w_mm, h_mm)

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
            ha='center', va='bottom',
            fontsize=style.body_pt() if p_value < 0.05 else style.tick_pt())

    # Styling
    ax.set_title(title)
    ax.set_ylabel(y_label)
    ax.set_xlabel('')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_ylim(ax.get_ylim()[0], bracket_top + y_range * 0.25)

    style.fit_margins(fig, pad_mm=0.6, cell_mm=cell,
                      reserve_letter=reserve_letter)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {w_mm} x {h_mm} mm canvas "
            f"(l,r,b,t mm): {over}")
    if reserve_letter:
        intruders = style.letter_clear(fig, cell)
        if intruders:
            raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    style.save_panel(fig, stem)
    print(f"  Saved: {stem.name}.[svg|pdf|png] at box {sub}, "
          f"{w_mm} x {h_mm} mm; letter cell {cell[0]:.2f} x {cell[1]:.2f} mm")
    print(f"    T: n={len(tumor_vals)}, mean={tumor_vals.mean():.4f}")
    print(f"    N: n={len(normal_vals)}, mean={normal_vals.mean():.4f}")
    print(f"    p-value: {p_value:.4e}")


def main():
    print("=" * 60)
    print("Panel H: External Validation Boxplots")
    print("=" * 60)

    script_dir = Path(__file__).parent

    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
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

    create_boxplot(df_linghua_filtered, 'GSE239676',
                   script_dir / 'mac3_boxplot_linghua', sub=1,
                   reserve_letter=True)

    # === Kumar Dataset ===
    print("\n[2/2] Kumar (GSE183904)...")
    df_kumar = pd.read_csv(FIG4_EXTERNAL_KUMAR)

    # Already has Normal vs Primary_tumor
    df_kumar['Group'] = df_kumar['Tissue'].map({
        'Primary_tumor': 'T',
        'Normal': 'N'
    })

    create_boxplot(df_kumar, 'GSE183904',
                   script_dir / 'mac3_boxplot_kumar', sub=2,
                   reserve_letter=False)

    print("\n" + "=" * 60)
    print("Panel H complete")
    print("=" * 60)


if __name__ == "__main__":
    main()
