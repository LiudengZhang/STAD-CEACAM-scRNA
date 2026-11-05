#!/usr/bin/env python3
"""Figure 4C: IM-MoMac module proportion in post-treatment stomach samples,
responders (Post-R) against non-responders (Post-NR). Two-sided Mann-Whitney
U test; the exact P is printed over the bracket.

Input:  FIG4_MODULE_PROPORTIONS, FIG4_CLINICAL_METADATA (00_Config/paths.py)
Output: module2_post_response_boxplot_k5.{svg,pdf,png} in this directory
Run:    python create_module2_boxplot.py
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
from scipy.stats import mannwhitneyu
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *                       # noqa: E402,F401,F403
from cnsfig import Panel, style  # noqa: E402
from cnsfig.boxes import draw_boxes, bracket, ylim_above, assert_no_points  # noqa: E402

PANEL = Panel(4, "C", __file__, small_pt=5.0)

# Response/phase colours shared by all figures.
from shared.figure_config import GROUP_COLORS, PRE_PAIR, POST_PAIR  # noqa: E402
PALETTE = list(POST_PAIR)   # post-treatment: the saturated pair, R first


def main():
    PANEL.apply()

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

    group_r = filtered[filtered['stomach_post_grouping'] == 'Responsed']['Module_2'].values
    group_nr = filtered[filtered['stomach_post_grouping'] == 'No-response']['Module_2'].values

    stat, p_value = mannwhitneyu(group_r, group_nr, alternative='two-sided')

    fig, ax = style.subplots_mm(PANEL.w_mm, PANEL.h_mm)

    bp = draw_boxes(ax, [group_r, group_nr], [0, 1], PALETTE, width=0.6)
    ax.set_xlim(-0.5, 1.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Post-R', 'Post-NR'])

    y_max = filtered['Module_2'].max()
    y_range = filtered['Module_2'].max() - filtered['Module_2'].min()
    bracket_top = y_max + y_range * 0.15
    _, p_txt, _ = bracket(fig, ax, 0, 1, y_max, y_range, p_value, kind="pair",
                          lift=0.10, arm=0.05, exact=True)

    ax.set_title('IM-MoMac')
    ax.set_ylabel('Module Proportion')
    ax.set_xlabel('')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_ylim(ax.get_ylim()[0], bracket_top + y_range * 0.15)
    ylim_above(ax, p_txt)
    assert_no_points(ax, bp)

    PANEL.fit(fig)
    PANEL.finish(fig, 'module2_post_response_boxplot_k5')


if __name__ == "__main__":
    main()
