#!/usr/bin/env python3
# REVISED FOR CIR-26-0753-ET, reviewer 1 point R1.3c:
# the test is two-sided and the annotation reports the exact P value.
"""
Figure 3 panel H - neighbourhood epithelial proportion in CEACAM-high against
CEACAM-low spots, paired within section.

The panel is drawn at the millimetre rectangle it prints in, read from
03_Final_Panels/panel_rects.csv through 00_Config/slots.py, so the type size
set here is the type size printed. Margins are measured from the rendered ink
rather than typed, and the panel-letter corner is left clear for the assembler.

  printed panel  Figure 3 H       (PROVENANCE.csv; the directory name agrees,
                                   but the letter was looked up, not inferred)

Three different printed panels - H, I and J - each write a file called
`ceacam_spatial_boxplot`. The stem is kept as it is: the assembler resolves
these by directory, and renaming one to disambiguate would break that.

Every value read, every filter and every statistic is the earlier drawing's.
The two-sided Wilcoxon signed-rank test, the sample-level means, the pairing,
the y limit and the bracket geometry are untouched.

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
    The published page prints "CEACAM-low" and "CEACAM-high" under the two
    boxes, on two lines. They were cut to "Low" and "High" on 2026-09-10
    because at 6 pt "CEACAM-" sets 9.76 mm while the two ticks stood between
    5.6 and 8.3 mm apart, so the pair overprinted.

    Restored 2026-09-11 in the two-line form the page uses. What changed is the
    room: these panels are 42.0 mm tall now rather than 31.0 (panel_rects_v2),
    a second label line costs height rather than width, and the axes no longer
    have to give up their own height to carry one. If the pair still overprints
    the panel gate says so - it is the check that convicted the single-line
    form - and the fallback is an axis label reading CEACAM region with Low and
    High on the ticks, not a silent second truncation.

    That is what happened. The two-line form was built and the gate convicted
    it on all three panels - 4.70, 0.48 and 1.99 mm2 of overlap, and on J a
    0.19 pt clearance against a y tick as well. The width is simply not there:
    these panels are 27.4 to 29.1 mm wide and that number did not change.

    So the ticks read Low and High and the axis is labelled CEACAM region. The
    qualifier is on the panel, under the boxes, where the published page puts
    it; it is set once instead of twice, which is the only reason it fits.

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
import panel_style_cns as style  # noqa: E402
import slots  # noqa: E402
from cnsfig.layout import pin_frame_mm  # noqa: E402
from cnsfig.boxes import box_xlim  # noqa: E402

SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 7.0                  # the earlier smallest body type, before * SCALE

MARK = style.tick_pt() / (SMALL_PT * SCALE)   # length multiplier
AREA = MARK ** 2                              # area multiplier

PANEL_LETTER = "H"
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
                    boxprops=dict(linewidth=style.RULE_PT),
                    whiskerprops=dict(linewidth=style.RULE_PT),
                    capprops=dict(linewidth=style.RULE_PT),
                    medianprops=dict(linewidth=style.RULE_PT),
                    flierprops=dict(markersize=6.0 * MARK,
                                    markeredgewidth=style.EDGE_PT))
    bp['boxes'][0].set_facecolor(CEACAM_LOW_COLOR)
    bp['boxes'][1].set_facecolor(CEACAM_HIGH_COLOR)
    for box in bp['boxes']:
        box.set_alpha(0.7)
        box.set_edgecolor('black')
    for median in bp['medians']:
        median.set_color('black')

    for i in range(len(low_data)):
        ax.plot([1, 2], [low_data[i], high_data[i]], 'k-', alpha=0.3, linewidth=style.RULE_PT)

    ax.scatter([1]*len(low_data), low_data, color=CEACAM_LOW_COLOR, s=40 * AREA, zorder=3,
               edgecolor='black', linewidth=style.EDGE_PT)
    ax.scatter([2]*len(high_data), high_data, color=CEACAM_HIGH_COLOR, s=40 * AREA, zorder=3,
               edgecolor='black', linewidth=style.EDGE_PT)

    data_y_max = max(max(low_data), max(high_data))
    y_range = data_y_max - min(min(low_data), min(high_data))
    bracket_y = data_y_max + 0.08 * y_range

    ax.plot([1, 1, 2, 2], [bracket_y, bracket_y + 0.02*y_range, bracket_y + 0.02*y_range, bracket_y],
            color='black', linewidth=style.RULE_PT)

    # A star at STAR_PT, a value at the tick size (panel_style_cns.p_text_kw,
    # 2026-09-15).
    pval_text, p_kw = style.p_text_kw(pval)
    ax.text(1.5, bracket_y + 0.04*y_range, pval_text, ha='center', va='bottom', **p_kw)

    # The title sits directly over the topmost tick label and at this size
    # their line boxes graze; a point of extra pad separates the two rows.
    ax.set_title(title, pad=plt.rcParams['axes.titlepad'] + 1.0)
    ax.set_ylabel(ylabel)
    ax.set_xticks([1, 2])
    # The page's two-line group names (2026-09-14). 'CEACAM-' sets 9.8 mm at
    # 6 pt, so the axis is widened at both ends to stand the ticks apart; the
    # gate measures the result.
    ax.set_xticklabels(['CEACAM-\nlow', 'CEACAM-\nhigh'])
    # Boxes off the spines (2026-09-14, evening): at xlim 0.75 the left box's
    # edge sat ON the y axis. 0.3 box widths of paper each side.
    ax.set_xlim(*box_xlim([1, 2], 0.5, clear=0.3))
    ax.tick_params(axis='both', width=style.RULE_PT, length=4 * MARK)
    # The lowest y tick label and the left group name met at the corner.
    ax.tick_params(axis='x', pad=2.5)
    # A shorter tick and a tighter pad on the y axis: the two group names
    # need every tenth of a millimetre this box can give them.
    ax.tick_params(axis='y', pad=0.6, length=1.5)

    for spine in ax.spines.values():
        spine.set_linewidth(style.RULE_PT)

    if y_max_limit is not None:
        ax.set_ylim(0, y_max_limit)

    # The sample count and the name of the test are stated in the caption; see
    # THE NOTE NAMING THE TEST above.
    print(f"  n = {len(low_data)} samples, Wilcoxon signed-rank (two-sided), "
          f"{pval_text}")

    # The title band above the plot clears the letter cell, so no left band
    # is reserved: the 3.3 mm goes to the plot, whose two group names need it.
    style.fit_margins(fig, pad_mm=0.3, cell_mm=LETTER_CELL, reserve_letter=False)
    # ONE FRAME LINE FOR F, G, H, I AND J  (2026-09-14, evening): the frame
    # top at 5.0 mm and its bottom at 29.0 mm below the slot top, so the five
    # plotting frames print 24 mm tall on one line (sweep_pages.ROW_ALIGN).
    pin_frame_mm(fig, ax, top_mm=5.0, bottom_mm=PANEL_H_MM - 29.0)
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
    print("Figure 3 panel H: neighbourhood epithelial proportion")
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
                    'neighborhood_epi_density': group_df['neighborhood_epi_density'].mean(),
                })

    sample_data = pd.DataFrame(sample_metrics)
    valid_samples = sample_data.groupby('sample').filter(lambda x: len(x) == 2)['sample'].unique()
    sample_data = sample_data[sample_data['sample'].isin(valid_samples)]
    sample_data = sample_data.sort_values(['sample', 'group'])
    print(f"  Samples with paired data: {len(valid_samples)}")

    print("\n[3/3] Creating boxplot...")
    create_paired_boxplot(sample_data, 'neighborhood_epi_density', 'Epithelial Density',
                         'Neighborhood\nEpithelial Proportion',
                         os.path.join(BASE_DIR, 'ceacam_spatial_boxplot'), y_max_limit=1.0)

    print("\nDone!")


if __name__ == '__main__':
    main()
