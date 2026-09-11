#!/usr/bin/env python3
# REVISED FOR CIR-26-0753-ET, reviewer 1 point R1.3c:
# the test is two-sided and the annotation reports the exact P value.
"""
Figure 3 panel J - distance to immune cells from CEACAM-high and CEACAM-low spots, paired
within section.

The panel is drawn at the millimetre rectangle it prints in, read from
03_Final_Panels/panel_rects.csv through 00_Config/slots.py, so the type size
set here is the type size printed. Margins are measured from the rendered ink
rather than typed, and the panel-letter corner is left clear for the assembler.

  printed panel  Figure 3 J       (PROVENANCE.csv; the directory name agrees,
                                   but the letter was looked up, not inferred)

Three different printed panels - H, I and J - each write a file called
`ceacam_spatial_boxplot`. The stem is kept as it is: the assembler resolves
these by directory, and renaming one to disambiguate would break that.

Every value read, every filter and every statistic is the earlier drawing's.
The two-sided Wilcoxon signed-rank test, the sample-level means, the pairing,
the y limit and the bracket geometry are untouched.

THE UNIT ON THE Y AXIS
    Until 2026-09-10 this panel plotted the raw `distance_to_distance_to_immune` column of
    `spot_data.csv` and labelled it `Distance (a.u.)`. That column is a
    Euclidean distance in the `x`, `y` of that table - full-resolution image
    pixels - and nothing in the path converted it, while the Results quoted the
    same quantity in micrometres. The figure was right and the sentence was not.

    Author's ruling 2026-09-10: convert both sides, so that the reader gets a
    physical distance and the page and the paragraph agree. The factor is the
    one the project already derives, `00_Config/spatial_scale.py` - the 100 um
    Visium centre-to-centre pitch over the median nearest-neighbour spacing of
    the array - and it is imported, never re-derived (RULES.md rule 5).

    ONE COHORT FACTOR, NOT TEN. The ten sections give 0.282476 to 0.288180
    um/unit, a 1.99 per cent spread. One constant for all ten makes this a pure
    linear rescale, so the paired Wilcoxon statistic and the exact P value
    printed on the panel are the ones that were printed before, to the digit;
    a per-section factor would re-weight the ten pairs against each other,
    which is an edit to the statistic and not to its unit. The reasoning is in
    `spatial_scale.cohort_um_per_unit`, which refuses to return one factor if
    the sections ever disagree by more than 5 per cent.

    The y limit is the old limit times the same factor, so the plotting box,
    the boxes, the pairing lines and the bracket sit exactly where they sat;
    only the numbers beside the axis and the word in its label change.

    THIS PANEL THEREFORE NO LONGER REPRODUCES THE PRINTED PAGE, deliberately
    and by ruling. See `SUPERSEDED.md` beside this script - the marker is not
    `KNOWN_BROKEN.md` and there is no `Retest-by` date, because there is
    nothing to retest.

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the corner note naming the test - at 7 * SCALE. MARK
    carries the non-type point sizes across to the 1:1 canvas:

        MARK = style.tick_pt() / (SMALL_PT * SCALE)

    The box, whisker, cap and median widths, the pairing lines, the dot area
    and its edge, the bracket and the tick width and length are scaled by it.
    Spine widths are not: those are style, and cnsplots sets them.

    The outlier marker is the one non-type point size the earlier drawing never
    named: it took matplotlib's default 6 pt marker and 1 pt edge on a canvas
    four times print size, so it printed at about 0.9 pt beside 4.3 pt type.
    cnsplots does not set it either, so it is rescaled by MARK like every other
    non-type length. The fliers themselves, and their values, are untouched.

THE GROUP TICK LABELS
    At 6 pt "CEACAM-" sets 9.76 mm while the two ticks stand between 5.6 and
    8.3 mm apart across the three panels, so the pair overprints; and the
    labels are themselves what holds the axes narrow, because the fit keeps
    them inside the canvas. The qualifier is the same on every box of all three
    panels and is stated in the caption, so each box keeps its own group name.
    Declared in RENAMES_FIGURE_3.

THE NOTE NAMING THE TEST
    The earlier drawing put "n = N samples / Wilcoxon signed-rank (two-sided)"
    inside the axes. At 6 pt the second line sets 31.6 mm and this panel is
    27.4 mm wide, so it cannot be drawn here at the figure's type size. The
    same sentence covers panels H, I and J, so it is stated once in the caption
    instead and declared in REMOVALS_FIGURE_3. The exact P value stays on the
    panel.
"""

import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import SPATIAL_SPOT_DATA
from spatial_scale import cohort_um_per_unit  # noqa: E402
import panel_style_cns as style  # noqa: E402
import slots  # noqa: E402

SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 7.0                  # the earlier smallest body type, before * SCALE

MARK = style.tick_pt() / (SMALL_PT * SCALE)   # length multiplier
AREA = MARK ** 2                              # area multiplier

PANEL_LETTER = "J"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(3, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(3, PANEL_LETTER)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

CEACAM_LOW_COLOR = '#1f77b4'
CEACAM_HIGH_COLOR = '#d62728'


def create_paired_boxplot(sample_data, col, title, ylabel, output_stem,
                          y_max_limit=None):
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
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

    # The title sits directly over the topmost tick label and at this size
    # their line boxes graze; a point of extra pad separates the two rows.
    ax.set_title(title, pad=plt.rcParams['axes.titlepad'] + 1.0)
    ax.set_ylabel(ylabel)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['Low', 'High'])
    ax.tick_params(axis='both', width=1.0 * MARK, length=4 * MARK)

    for spine in ax.spines.values():
        spine.set_linewidth(1.0 * MARK)

    if y_max_limit is not None:
        ax.set_ylim(0, y_max_limit)

    # The sample count and the name of the test are stated in the caption; see
    # THE NOTE NAMING THE TEST above.
    print(f"  n = {len(low_data)} samples, Wilcoxon signed-rank (two-sided), "
          f"{pval_text}")

    style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_CELL)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    style.save_panel(fig, output_stem)
    print(f"  Saved: {output_stem}.[svg|pdf|png] at "
          f"{PANEL_W_MM} x {PANEL_H_MM} mm")


def main():
    print("=" * 60)
    print("Figure 3 panel J: distance to immune")
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

    # Array units -> micrometres, one factor for the cohort, derived in
    # 00_Config/spatial_scale.py and imported rather than repeated. See
    # THE UNIT ON THE Y AXIS above. A single constant is a rescale: the
    # paired Wilcoxon statistic and the printed P value cannot move.
    um_per_unit = cohort_um_per_unit(df)
    print(f"  distances -> micrometres at {um_per_unit:.9f} um per array unit")
    sample_data['distance_to_immune'] = sample_data['distance_to_immune'] * um_per_unit
    print(f"  Samples with paired data: {len(valid_samples)}")

    print("\n[3/3] Creating boxplot...")
    create_paired_boxplot(sample_data, 'distance_to_immune', 'Distance to Immune',
                         'Distance (µm)',
                         os.path.join(BASE_DIR, 'ceacam_spatial_boxplot'),
                         y_max_limit=1200 * um_per_unit)

    print("\nDone!")


if __name__ == '__main__':
    main()
