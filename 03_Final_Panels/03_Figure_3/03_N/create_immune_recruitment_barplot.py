#!/usr/bin/env python3
"""
Figure 3N: difference in immune-cell proportion between CEACAM5/6-high and
CEACAM5/6-low regions, one bar per cell type (red above zero, blue below),
with stars for the per-cell-type P values in the table (uncorrected).

Reads SPATIAL_REGION_COMPARISON (00_Config/paths.py).
Writes ceacam_immune_region_summary.[svg|pdf|png] beside this script.

Run: python create_immune_recruitment_barplot.py
"""

import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import SPATIAL_REGION_COMPARISON
from cnsfig import Panel, pin_frame_mm, style  # noqa: E402

PANEL = Panel(3, "N", __file__, small_pt=7.0)

# The y label on three lines; on two it would be taller than the frame.
YLABEL = 'Difference in Proportion\n(CEACAM5/6-high -\nCEACAM5/6-low)'
# Gap between a bar's end and its star's ink.
STAR_GAP_MM = 0.4
MARK, AREA = PANEL.mark, PANEL.area


def main():
    results_df = pd.read_csv(SPATIAL_REGION_COMPARISON)
    results_df = results_df.sort_values('difference', ascending=False)

    fig, ax = PANEL.subplots()

    SHORT_NAMES = {
        'Monocytes/Macrophages': 'Mo/Mac',
        'CD4+ T cells': 'CD4+ T',
        'CD8+ T cells': 'CD8+ T',
        'Plasma cells': 'Plasma',
    }

    df_sorted = results_df.sort_values('difference', ascending=False)
    df_sorted['label'] = df_sorted['immune_cell'].map(lambda x: SHORT_NAMES.get(x, x))
    x_pos = np.arange(len(df_sorted))
    colors = ['#d62728' if d > 0 else '#1f77b4' for d in df_sorted['difference']]

    bars = ax.bar(x_pos, df_sorted['difference'], color=colors, edgecolor='black',
                  alpha=0.8, width=0.7, linewidth=style.RULE_PT)

    ax.axhline(y=0, color='black', linestyle='-', linewidth=style.RULE_PT)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(df_sorted['label'], rotation=45, ha='right',
                       rotation_mode='anchor')
    ax.set_ylabel(YLABEL)
    ax.set_title('Immune Composition by CEACAM5/6 Region', fontweight='normal')
    ax.tick_params(axis='both', width=style.RULE_PT, length=4 * MARK)
    for spine in ax.spines.values():
        spine.set_linewidth(style.RULE_PT)

    # Each star sits STAR_GAP_MM beyond its bar's end, measured from the
    # glyph's ink extent.
    from matplotlib import transforms
    from matplotlib.textpath import TextPath
    from matplotlib.font_manager import FontProperties
    gap_pt = STAR_GAP_MM * 72 / 25.4
    y_min = df_sorted['difference'].min() - 0.022
    y_max = df_sorted['difference'].max() + 0.02
    ax.set_ylim(y_min, y_max)
    for i, (_, row) in enumerate(df_sorted.iterrows()):
        pval = row['pvalue']
        if not style.p_is_star(pval):
            continue
        sig, kw = style.p_text_kw(pval)
        ext = TextPath((0, 0), sig, size=kw['fontsize'],
                       prop=FontProperties(family=style.resolved_family())
                       ).get_extents()
        diff = row['difference']
        dy = (gap_pt - ext.y0) if diff >= 0 else (-gap_pt - ext.y1)
        tr = transforms.offset_copy(ax.transData, fig=fig, x=0, y=dy,
                                    units='points')
        ax.text(i, diff, sig, ha='center', va='baseline', transform=tr, **kw)

    PANEL.fit(fig)
    # Same frame top and height as the spatial maps in this row.
    pin_frame_mm(fig, ax, top_mm=5.0, bottom_mm=PANEL.h_mm - 29.0)
    PANEL.finish(fig, 'ceacam_immune_region_summary')


if __name__ == '__main__':
    main()
