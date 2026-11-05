#!/usr/bin/env python3
"""
Figure 2F: proportion of the Epi_CEACAM5/6 cluster (C2) among pre-treatment
epithelial cells per sample, responders versus non-responders.

Reads the epithelial h5ad (EPITHELIAL_H5AD).
Writes c2_proportion_pre_boxplot.[svg|pdf|png] beside this script.

Run: python create_c2_proportion_boxplot.py

The panel prints the exact one-tailed Mann-Whitney P (non-responders >
responders); the two-sided P is reported in the legend.
"""

import scanpy as sc
import numpy as np
import os

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import EPITHELIAL_H5AD                          # noqa: E402
from cnsfig import Panel, style  # noqa: E402
from cnsfig.boxes import draw_boxes, bracket, ylim_above, assert_no_points, printed_and_two_sided  # noqa: E402

# Shared response/phase colours.
from shared.figure_config import GROUP_COLORS, PRE_PAIR, POST_PAIR  # noqa: E402
# Pre-treatment samples: the light pair.
COLORS = {'Responsed': PRE_PAIR[0], 'No-response': PRE_PAIR[1]}

SCALE = 4                           # sets MARK relative to the tick size
SMALL_PT = 6.0

PANEL = Panel(2, "F", __file__, small_pt=SMALL_PT)


def main():
    style.apply(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)

    adata = sc.read_h5ad(EPITHELIAL_H5AD)

    pre_mask = adata.obs['Treatment phase'] == 'Pre'
    valid_mask = adata.obs['stomach_pre_grouping'].isin(['Responsed', 'No-response'])
    adata_filtered = adata[pre_mask & valid_mask].copy()

    adata_filtered.obs['is_C2'] = (adata_filtered.obs['minor_cell_state'] == 'C2_Epi_CEACAM6').astype(int)

    sample_props = adata_filtered.obs.groupby('sample', observed=True).agg({
        'is_C2': 'mean',
        'stomach_pre_grouping': 'first'
    }).reset_index()
    sample_props['C2_proportion'] = sample_props['is_C2'] * 100

    responder_vals = sample_props[sample_props['stomach_pre_grouping'] == 'Responsed']['C2_proportion'].values
    non_responder_vals = sample_props[sample_props['stomach_pre_grouping'] == 'No-response']['C2_proportion'].values

    # One-tailed P for the panel, two-sided P for the legend, from one call.
    pval, pval_two_sided = printed_and_two_sided(non_responder_vals, responder_vals)

    fig, ax = style.subplots_mm(PANEL.w_mm, PANEL.h_mm)

    data = [responder_vals, non_responder_vals]
    bp = draw_boxes(ax, data, [1, 2],
                    [COLORS['Responsed'], COLORS['No-response']], width=0.6)

    y_max = max(np.max(responder_vals), np.max(non_responder_vals))
    _, p_text, _ = bracket(fig, ax, 1, 2, y_max, y_max, pval, kind="pair",
                           lift=0.15, arm=0.05 * 1.15, exact=True)

    ax.set_title('Epi_CEACAM5/6')
    ax.set_ylabel('Proportion (%)')
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['Pre-R', 'Pre-NR'])
    ax.set_ylim(0, y_max * 1.40)
    ylim_above(ax, p_text)
    assert_no_points(ax, bp)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    PANEL.fit(fig)

    os.makedirs(PANEL.dir, exist_ok=True)
    PANEL.finish(fig, 'c2_proportion_pre_boxplot')


if __name__ == '__main__':
    main()
