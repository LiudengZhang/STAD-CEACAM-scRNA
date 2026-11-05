#!/usr/bin/env python3
"""
Figure 2K, left box: sample-mean CEACAM6 expression in pre-treatment stomach
epithelium, responders versus non-responders. draw() is shared with the right
box (02_K/create_ceacam5_boxplot.py, CEACAM5).

The panel prints the exact one-tailed Mann-Whitney P (non-responders >
responders); the two-sided P is reported in the legend.

Reads the epithelial h5ad (EPITHELIAL_H5AD); expression from .raw, or .X when
the file has no .raw.
Writes ceacam6_pre_boxplot.[svg|pdf|png] beside this script.

Run: python create_ceacam6_boxplot.py
"""

import os
import sys
from pathlib import Path

from matplotlib.ticker import MultipleLocator
import numpy as np
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import EPITHELIAL_H5AD                                    # noqa: E402
from shared.expression import assert_x_stands_in_for_raw             # noqa: E402
from shared.figure_config import RESPONSE_COLORS                     # noqa: E402
import panel_style_cns as style                                      # noqa: E402
import slots                                                         # noqa: E402
from cnsfig.boxes import box_xlim, frame_fixed, draw_boxes, bracket, ylim_above, assert_no_points, printed_and_two_sided  # noqa: E402

# Shared response/phase colours.
from shared.figure_config import GROUP_COLORS, PRE_PAIR, POST_PAIR  # noqa: E402
# Pre-treatment samples: the light pair.
COLORS = {'Responsed': PRE_PAIR[0], 'No-response': PRE_PAIR[1]}

PANEL_LETTER = "K"

SCALE = 4
SMALL_PT = 6.0


def draw(out_dir, sub, gene, stem, ylabel, left_mm):
    """One gene's R-versus-NR box in sub-box `sub` of Figure 2K.

    The four boxes of Figure 2K and 2L share fixed mm margins (frame_fixed;
    left 9.0 mm with a y label, 5.4 mm without), so their frames line up.
    printed_and_two_sided returns the one-tailed P the panel prints and the
    two-sided P for the legend from one call.
    """
    panel_w_mm, panel_h_mm = slots.size_mm(2, PANEL_LETTER, sub=sub)

    style.apply(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)

    adata = sc.read_h5ad(EPITHELIAL_H5AD)
    # Raises if this file has no .raw and its .X is not known to equal it.
    assert_x_stands_in_for_raw(adata, EPITHELIAL_H5AD)

    pre_mask = adata.obs['Treatment phase'] == 'Pre'
    valid_mask = adata.obs['stomach_pre_grouping'].isin(['Responsed', 'No-response'])
    adata_filtered = adata[pre_mask & valid_mask].copy()

    # Log1p expression from .raw; the deposited h5ad has no .raw and holds the
    # same matrix in .X.
    src = adata_filtered.raw if adata_filtered.raw is not None else adata_filtered
    idx = src.var_names.get_loc(gene)
    adata_filtered.obs[gene] = src.X[:, idx].toarray().flatten()

    sample_means = adata_filtered.obs.groupby('sample', observed=True).agg({
        gene: 'mean',
        'stomach_pre_grouping': 'first'
    }).reset_index()

    responder_vals = sample_means[sample_means['stomach_pre_grouping'] == 'Responsed'][gene].values
    non_responder_vals = sample_means[sample_means['stomach_pre_grouping'] == 'No-response'][gene].values

    pval, pval_two_sided = printed_and_two_sided(non_responder_vals, responder_vals)

    fig, ax = style.subplots_mm(panel_w_mm, panel_h_mm)

    data = [responder_vals, non_responder_vals]
    positions = [1, 2]
    bp = draw_boxes(ax, data, positions,
                    [COLORS['Responsed'], COLORS['No-response']], width=0.6)
    y_max = max(np.max(responder_vals), np.max(non_responder_vals))
    # Bracket at 1.15 x y_max with the exact P above it.
    _, p_text, _ = bracket(fig, ax, 1, 2, y_max, y_max, pval, kind="pair",
                           lift=0.15, arm=0.05 * 1.15, exact=True)

    ax.set_xticks(positions)
    ax.set_xticklabels(['Pre-R', 'Pre-NR'])
    ax.set_xlim(*box_xlim(positions, 0.6, clear=0.3))   # narrow frame: room for the group names
    ax.set_ylim(0, y_max * 1.75)   # headroom: the P string clears the title
    ylim_above(ax, p_text)
    assert_no_points(ax, bp)
    # y ticks every 0.5.
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    frame_fixed(fig, ax, title=f'*{gene}*', ylabel=ylabel,
                left_mm=left_mm, right_mm=0.6, bottom_mm=3.8,
                panel_w_mm=panel_w_mm, panel_h_mm=panel_h_mm)

    os.makedirs(out_dir, exist_ok=True)
    return style.save_panel(fig, Path(out_dir) / stem)


def main():
    draw(Path(__file__).parent, 1, 'CEACAM6', 'ceacam6_pre_boxplot',
         'Expression', 9.0)


if __name__ == '__main__':
    main()
