#!/usr/bin/env python3
# REVISED FOR CIR-26-0753-ET, reviewer 1 point R1.3c:
# the test is two-sided and the annotation reports the exact P value.
# The one-tailed version this replaces is the one used in the preprint,
# https://www.biorxiv.org/content/10.64898/2026.03.05.708917
"""
Figure 2 panel K, CEACAM6 box - sample-mean CEACAM6 expression in
pre-treatment stomach epithelium, responders versus non-responders.

Printed K is one letter over two drawings, 02_J (CEACAM6, this file) and 02_K (CEACAM5). Each is
drawn at its own printed sub-box, read from 03_Final_Panels/slot_subrects.csv
through 00_Config/slots.py, so the type size set here is the type size printed.
Margins are measured from the rendered ink rather than typed.

The printed letter K sits to the left of both sub-boxes, so neither of them
reserves a corner for it.

  printed panel  Figure 2 K, box 1   (PROVENANCE.csv; NOT inferred from "02_J")

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the y axis label, the tick labels and the bracket P
    value - at 6 * SCALE. MARK carries the non-type point sizes across to the
    1:1 canvas:

        MARK = style.tick_pt() / (SMALL_PT * SCALE)

    Every non-type length - the box, whisker, cap, median and bracket line
    widths and the marker sizes - takes it. Tick widths and lengths and spine
    widths do not: those are style, and cnsplots sets them.

Every value read, every filter, every statistic and every string is the earlier
drawing's. The drawing code is the same code.
"""

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import scanpy as sc
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import EPITHELIAL_H5AD                                    # noqa: E402
from shared.figure_config import (RESPONSE_COLORS,                   # noqa: E402
                                  RESPONSE_MEDIAN_COLORS)
import panel_style_cns as style                                      # noqa: E402
import slots                                                         # noqa: E402
from cnsfig.boxes import box_xlim, frame_fixed, draw_boxes, bracket, ylim_above, assert_no_points  # noqa: E402

COLORS = RESPONSE_COLORS
MEDIAN_COLORS = RESPONSE_MEDIAN_COLORS

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = BASE_DIR

PANEL_LETTER = "K"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(2, PANEL_LETTER, sub=1)

SCALE = 4                           # the earlier canvas multiplier
SMALL_PT = 6.0                      # the earlier smallest body type


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.5f}")

    # Load data
    print("Loading data...")
    adata = sc.read_h5ad(EPITHELIAL_H5AD)
    print(f"Loaded {adata.n_obs} cells")

    # Filter for PRE-TREATMENT only with valid response groups
    pre_mask = adata.obs['Treatment phase'] == 'Pre'
    valid_mask = adata.obs['stomach_pre_grouping'].isin(['Responsed', 'No-response'])
    adata_filtered = adata[pre_mask & valid_mask].copy()
    print(f"Pre-treatment cells with response status: {adata_filtered.n_obs}")

    # Expression comes from .raw wherever there is a .raw: .X in the working
    # inputs is the scaled matrix left by a double normalisation and carries
    # NaN rows. The clean deposit promotes .raw.X to .X and carries no .raw at
    # all, so a file without one holds the same log1p numbers in .X.
    src = adata_filtered.raw if adata_filtered.raw is not None else adata_filtered
    ceacam6_idx = src.var_names.get_loc('CEACAM6')
    ceacam6_expr = src.X[:, ceacam6_idx].toarray().flatten()

    # Add to obs for aggregation
    adata_filtered.obs['CEACAM6'] = ceacam6_expr

    # Calculate sample-level mean expression
    sample_means = adata_filtered.obs.groupby('sample', observed=True).agg({
        'CEACAM6': 'mean',
        'stomach_pre_grouping': 'first'
    }).reset_index()
    print(f"\nSample-level data ({len(sample_means)} samples):")

    # Separate by response group
    responder_vals = sample_means[sample_means['stomach_pre_grouping'] == 'Responsed']['CEACAM6'].values
    non_responder_vals = sample_means[sample_means['stomach_pre_grouping'] == 'No-response']['CEACAM6'].values

    print(f"\nResponders (n={len(responder_vals)}): mean={np.mean(responder_vals):.3f}")
    print(f"Non-Responders (n={len(non_responder_vals)}): mean={np.mean(non_responder_vals):.3f}")

    # Mann-Whitney U test
    stat, pval = stats.mannwhitneyu(non_responder_vals, responder_vals, alternative='two-sided')
    print(f"\nTwo-sided Mann-Whitney U test: U={stat:.2f}, p={pval:.4f}")

    # Create figure at the millimetre size it prints at
    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

    # Prepare data for boxplot
    data = [responder_vals, non_responder_vals]
    positions = [1, 2]

    # ONE BOX (2026-09-16, the author's fifth reading: "box plots should
    # look alike throughout"): cnsfig.boxes.draw_boxes - the 2H/I box, black
    # median (the per-group median colour is the one thing given up). The
    # bracket keeps its vertices (line at 1.15 x y_max, arms 5% of that);
    # the P string's ink sits 0.4 mm above the line (cnsfig.boxes.bracket).
    bp = draw_boxes(ax, data, positions,
                    [COLORS['Responsed'], COLORS['No-response']], width=0.6)
    y_max = max(np.max(responder_vals), np.max(non_responder_vals))
    _, p_text, _ = bracket(fig, ax, 1, 2, y_max, y_max, pval, kind="pair",
                           lift=0.15, arm=0.05 * 1.15)

    # Labels - sizes now come from cnsplots' rcParams
    # The gene symbol is set in italic, as the shipped page sets it and as
    # the rest of this figure sets one. Style only; the string is unchanged.
    # ONE FRAME FOR THE FOUR BOXES OF K AND L  (2026-09-14, evening)
    #   The author asked for K and L to align with each other and with J.
    #   The four boxes name the same margins in millimetres (cnsfig.boxes
    #   frame_fixed: left 9.0 mm for a box with a y label, 5.4 mm without;
    #   3.8 mm below for the tick labels; one title line above), so their
    #   frames print at the same x, the pairs at the same y, and the x
    #   limits come from box_xlim so the boxes stand off the spines.
    ax.set_xticks(positions)
    ax.set_xticklabels(['Pre-R', 'Pre-NR'])
    ax.set_xlim(*box_xlim(positions, 0.6, clear=0.3))   # 14 mm frame: the two group names need the pitch
    ax.set_ylim(0, y_max * 1.75)   # headroom: the P string clears the title
    ylim_above(ax, p_text)
    assert_no_points(ax, bp)
    # The published ruler: 0.0 to 2.0 in steps of 0.5 (2026-09-14).
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    frame_fixed(fig, ax, title='*CEACAM6*', ylabel='Expression',
                left_mm=9.0, right_mm=0.6, bottom_mm=3.8,
                panel_w_mm=PANEL_W_MM, panel_h_mm=PANEL_H_MM)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    style.save_panel(fig, Path(OUTPUT_DIR) / "ceacam6_pre_boxplot")
    print(f"\nSaved: {Path(OUTPUT_DIR) / 'ceacam6_pre_boxplot'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
