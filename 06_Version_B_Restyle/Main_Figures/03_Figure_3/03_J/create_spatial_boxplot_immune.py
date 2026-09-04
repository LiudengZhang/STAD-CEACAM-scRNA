#!/usr/bin/env python3
# REVISED FOR CIR-26-0753-ET, reviewer 1 point R1.3c:
# the test is two-sided and the annotation reports the exact P value.
# Original one-tailed version: Round_5/03_Final_Panels/03_Figure_3/03_J/create_spatial_boxplot_immune.py
"""
Figure 3 panel J, RESTYLED (Version B) - distance to immune, CEACAM-high vs
CEACAM-low, paired within section.

NOTE THE PANEL LETTER, AND THE FILENAME. This directory is `03_J` and
PROVENANCE.csv says it holds printed panel **J** (CLAUDE.md rule 2: the
letter is looked up, never inferred). Three *different* printed panels - H, I
and J - each write a file called `ceacam_spatial_boxplot.png`, in `03_H`,
`03_I` and `03_J` respectively. The stem is kept exactly as Version A has it:
the assembler resolves these by directory, and renaming one to disambiguate
would break that.

Version A is
`03_Final_Panels/03_Figure_3/03_J/create_spatial_boxplot_immune.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every value read, every filter and every statistic is Version A's. The
two-sided Wilcoxon signed-rank test, the sample-level means, the pairing, the
`y_max_limit=1200` and the bracket geometry are untouched.

  printed panel  Figure 3 J     (PROVENANCE.csv; NOT inferred from "03_J")
  printed rect   28.6 x 32.1 mm   (panel_rects.csv)
  Version B box  52.0 x 54.0 mm

MARK
    Version A drew a 20.0 x 22.0 cm canvas (SCALE = 4). Its smallest body type
    is the corner note - "n = N samples / Wilcoxon signed-rank (two-sided)" -
    at `fontsize=7 * SCALE`, so SMALL_PT = 7 and, by PANEL_SPEC.md,

        MARK = tick_pt / (SMALL_PT * SCALE) = 7 / 28 = 0.25
        AREA = MARK ** 2                    = 0.0625

    Every non-type length Version A set was a literal point size on that 4x
    canvas - box/whisker/cap 1.0, median 2.0, pairing lines 0.8, dot area 40
    with a 0.5 edge, the bracket 1, the spines 1.0 and the ticks 1.0 / 4 - and
    each is multiplied by it, so all keep their Version A size *relative to the
    type*.

    That corner note is also why the box is 52.0 mm wide and not 28.6: at 7 pt
    "Wilcoxon signed-rank (two-sided)" sets ~28 mm, and Version A placed it
    inside the axes at `transAxes` x = 0.98. The axes has to be wider than the
    note or the note runs out of the panel. Nothing plotted moved.
"""

import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import SPATIAL_SPOT_DATA
import panel_style_cns as style  # noqa: E402

# Version A's canvas multiplier. Version B draws 1:1, so SCALE survives only to
# reproduce the exact numbers Version A set for its non-type lengths.
SCALE = 4
SMALL_PT = 7.0                       # Version A's smallest body type

# See MARK above. `style.tick_pt()` reads cnsplots' own setting, so the factor
# is derived, never a literal.
MARK = style.tick_pt() / (SMALL_PT * SCALE)   # 0.25, length multiplier
AREA = MARK ** 2                              # 0.0625, area multiplier

PRINTED_MM = (28.6, 32.1)          # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 52.0, 54.0
# Millimetres of paper. left holds the y label and its tick labels, bottom the
# two-line group labels, top the title.
MARGIN = dict(left=13.0, right=2.0, top=5.0, bottom=9.0)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

CEACAM_LOW_COLOR = '#1f77b4'
CEACAM_HIGH_COLOR = '#d62728'


def create_paired_boxplot(sample_data, col, title, ylabel, output_stem, y_max_limit=None):
    family = style.apply()
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

    low_data = sample_data[sample_data['group'] == 'CEACAM-low'][col].values
    high_data = sample_data[sample_data['group'] == 'CEACAM-high'][col].values

    stat, pval = stats.wilcoxon(high_data, low_data, alternative='two-sided')

    bp = ax.boxplot([low_data, high_data], positions=[1, 2], widths=0.5, patch_artist=True,
                    boxprops=dict(linewidth=1.0 * MARK),
                    whiskerprops=dict(linewidth=1.0 * MARK),
                    capprops=dict(linewidth=1.0 * MARK),
                    medianprops=dict(linewidth=2.0 * MARK),
                    # The outlier marker is the one non-type point size Version A
                    # never named: it took matplotlib's rcParams default of 6 pt
                    # markersize / 1 pt edge on a canvas four times print size, so
                    # it printed at ~0.9 pt beside ~4.3 pt type. cnsplots does not
                    # set it either, so leaving it alone would print it at 6 pt
                    # beside 7 pt type - six times too big relative to everything
                    # around it. It is rescaled by MARK like every other non-type
                    # length. The fliers themselves, and their values, are
                    # untouched.
                    flierprops=dict(markersize=6.0 * MARK,
                                    markeredgewidth=1.0 * MARK))
    bp['boxes'][0].set_facecolor(CEACAM_LOW_COLOR)
    bp['boxes'][1].set_facecolor(CEACAM_HIGH_COLOR)
    for box in bp['boxes']:
        box.set_alpha(0.7)
        box.set_edgecolor('black')
    for median in bp['medians']:
        median.set_color('black')

    for i in range(len(low_data)):
        ax.plot([1, 2], [low_data[i], high_data[i]], 'k-', alpha=0.3, linewidth=0.8 * MARK)

    ax.scatter([1]*len(low_data), low_data, color=CEACAM_LOW_COLOR, s=40 * AREA, zorder=3,
               edgecolor='black', linewidth=0.5 * MARK)
    ax.scatter([2]*len(high_data), high_data, color=CEACAM_HIGH_COLOR, s=40 * AREA, zorder=3,
               edgecolor='black', linewidth=0.5 * MARK)

    data_y_max = max(max(low_data), max(high_data))
    y_range = data_y_max - min(min(low_data), min(high_data))
    bracket_y = data_y_max + 0.08 * y_range

    ax.plot([1, 1, 2, 2], [bracket_y, bracket_y + 0.02*y_range, bracket_y + 0.02*y_range, bracket_y],
            color='black', linewidth=1 * MARK)

    pval_text = f'P = {pval:.3f}' if pval >= 0.001 else 'P < 0.001'
    ax.text(1.5, bracket_y + 0.04*y_range, pval_text, ha='center', va='bottom')

    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['CEACAM-\nlow', 'CEACAM-\nhigh'])
    ax.tick_params(axis='both', width=1.0 * MARK, length=4 * MARK)

    for spine in ax.spines.values():
        spine.set_linewidth(1.0 * MARK)

    if y_max_limit is not None:
        ax.set_ylim(0, y_max_limit)

    n_pairs = len(low_data)
    ax.text(0.98, 0.02, f'n = {n_pairs} samples\nWilcoxon signed-rank (two-sided)',
            transform=ax.transAxes, fontsize=style.tick_pt(), va='bottom', ha='right',
            color='#555555')

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, output_stem)

    print(f"  Saved: {output_stem}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")


def main():
    print("=" * 60)
    print("Panel G: Distance to Immune Boxplot")
    print("=" * 60)

    print("\n[1/3] Loading spot data...")
    df = pd.read_csv(SPATIAL_SPOT_DATA)
    print(f"  Total spots (all samples): {len(df):,}")

    print("\n[2/3] Computing sample-level means...")
    sample_metrics = []
    for sample in df['sample'].unique():
        sample_df = df[df['sample'] == sample]
        for group in ['CEACAM-high', 'CEACAM-low']:
            group_df = sample_df[sample_df['CEACAM_group'] == group]
            if len(group_df) > 0:
                sample_metrics.append({
                    'sample': sample,
                    'group': group,
                    'distance_to_immune': group_df['distance_to_immune'].mean(),
                })

    sample_data = pd.DataFrame(sample_metrics)
    valid_samples = sample_data.groupby('sample').filter(lambda x: len(x) == 2)['sample'].unique()
    sample_data = sample_data[sample_data['sample'].isin(valid_samples)]
    sample_data = sample_data.sort_values(['sample', 'group'])
    print(f"  Samples with paired data: {len(valid_samples)}")

    print("\n[3/3] Creating boxplot...")
    create_paired_boxplot(sample_data, 'distance_to_immune', 'Distance to Immune',
                         'Distance (a.u.)',
                         os.path.join(BASE_DIR, 'ceacam_spatial_boxplot'), y_max_limit=1200)

    print("\nDone!")


if __name__ == '__main__':
    main()
