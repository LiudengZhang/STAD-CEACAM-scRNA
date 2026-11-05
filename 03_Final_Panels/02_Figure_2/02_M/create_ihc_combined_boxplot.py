#!/usr/bin/env python3
"""
Figure 2N: combined CEACAM5 + CEACAM6 immunohistochemical staining area in
the eight pre-treatment specimens, responders versus non-responders.

Staining area comes from colour deconvolution (Ruifrok & Johnston 2001;
skimage rgb2hed, DAB optical-density threshold 0.02), computed by
upstream/IHC/quantify_ceacam_ihc.py. The panel prints the
exact one-tailed Mann-Whitney P (non-responders > responders); the two-sided
P is reported in the legend. Individual points are not drawn.

Reads IHC_COLOR_DECONV_CSV (00_Config/paths.py).
Writes ihc_combined_boxplot.[svg|pdf|png] beside this script.

Run: python create_ihc_combined_boxplot.py
"""

import numpy as np
import pandas as pd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))

# Shared response/phase colours, as in Figure 2K.
from shared.figure_config import GROUP_COLORS, PRE_PAIR, POST_PAIR  # noqa: E402
COLOR_R, COLOR_NR = PRE_PAIR          # pre-treatment IHC: the light pair

SCALE = 4                           # sets MARK relative to the tick size
SMALL_PT = 6.0

from paths import IHC_COLOR_DECONV_CSV  # noqa: E402
from cnsfig import Panel, style  # noqa: E402
from cnsfig.boxes import draw_boxes, bracket, ylim_above, assert_no_points, printed_and_two_sided  # noqa: E402

PANEL = Panel(2, "N", __file__, small_pt=SMALL_PT)
IHC_CSV = IHC_COLOR_DECONV_CSV
df = pd.read_csv(IHC_CSV)

# Pivot: combined CEACAM5 + CEACAM6 per patient
pivot = df.pivot_table(index=['patient', 'group'], columns='marker',
                       values='staining_pct').reset_index()
pivot['combined'] = pivot['CEACAM5'] + pivot['CEACAM6']

R_VALS = pivot[pivot['group'] == 'R']['combined'].values
NR_VALS = pivot[pivot['group'] == 'NR']['combined'].values


def main():
    style.apply(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2

    fig, ax = style.subplots_mm(PANEL.w_mm, PANEL.h_mm)

    # One-tailed P for the panel, two-sided P for the legend, from one call.
    pval, pval_two_sided = printed_and_two_sided(NR_VALS, R_VALS)

    bp = draw_boxes(ax, [R_VALS, NR_VALS], [1, 2], [COLOR_R, COLOR_NR],
                    width=0.6)

    y_max = max(np.max(R_VALS), np.max(NR_VALS))
    # Bracket at 1.15 x y_max with the exact P above it.
    _, p_text, _ = bracket(fig, ax, 1, 2, y_max, y_max, pval, kind="pair",
                           lift=0.15, arm=0.05 * 1.15, exact=True)

    ax.set_title('CEACAM5/6 IHC')
    ax.set_ylabel('Staining area (%)')
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['Pre-R', 'Pre-NR'])
    ax.set_ylim(0, y_max * 1.40)
    ylim_above(ax, p_text)
    assert_no_points(ax, bp)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    PANEL.fit(fig)

    PANEL.finish(fig, 'ihc_combined_boxplot')


if __name__ == '__main__':
    main()
