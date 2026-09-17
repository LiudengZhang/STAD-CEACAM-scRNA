#!/usr/bin/env python3
"""
Figure 2 panel E - C2_Epi_CEACAM5/6 cluster proportion in pre-treatment
stomach epithelium, responders versus non-responders.

The panel is drawn at the millimetre rectangle it prints in, read from
03_Final_Panels/panel_rects.csv through 00_Config/slots.py, so the type size
set here is the type size printed. Margins are measured from the rendered ink
rather than typed, and the panel-letter corner is left clear for the assembler.

  printed panel  Figure 2 E       (PROVENANCE.csv; the directory name agrees,
                                   but the letter was looked up, not inferred)

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the y axis label, the tick labels and the bracket P
    value - at 6 * SCALE. MARK carries the non-type point sizes across to the
    1:1 canvas:

        MARK = style.tick_pt() / (SMALL_PT * SCALE)

    The box, whisker, cap, median and bracket line widths and the flier marker
    take it. Tick widths and lengths and spine widths do not: those are style,
    and cnsplots sets them.

Every value read, every filter, every statistic and every string is the earlier
drawing's. The drawing code is the same code.
"""

import scanpy as sc
import numpy as np
from scipy import stats
import os

# Central config
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *                                       # noqa: E402,F403
import panel_style_cns as style                           # noqa: E402
import slots                                              # noqa: E402
from cnsfig.boxes import draw_boxes, bracket, ylim_above, assert_no_points  # noqa: E402

# Paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = BASE_DIR

# Colors: Blue for R, Red for NR
COLORS = {
    'Responsed': '#0072B2',
    'No-response': '#D55E00',
}
MEDIAN_COLORS = {
    'Responsed': '#005689',
    'No-response': '#A34700',
}

PANEL_LETTER = "E"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(2, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(2, PANEL_LETTER)

SCALE = 4                           # the earlier canvas multiplier, for MARK
SMALL_PT = 6.0                      # the earlier smallest body type


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.5f}")

    print("Loading data...")
    adata = sc.read_h5ad(EPITHELIAL_H5AD)                 # noqa: F405

    pre_mask = adata.obs['Treatment phase'] == 'Pre'
    valid_mask = adata.obs['stomach_pre_grouping'].isin(['Responsed', 'No-response'])
    adata_filtered = adata[pre_mask & valid_mask].copy()
    print(f"Pre-treatment cells: {adata_filtered.n_obs}")

    adata_filtered.obs['is_C2'] = (adata_filtered.obs['minor_cell_state'] == 'C2_Epi_CEACAM6').astype(int)

    sample_props = adata_filtered.obs.groupby('sample', observed=True).agg({
        'is_C2': 'mean',
        'stomach_pre_grouping': 'first'
    }).reset_index()
    sample_props['C2_proportion'] = sample_props['is_C2'] * 100

    responder_vals = sample_props[sample_props['stomach_pre_grouping'] == 'Responsed']['C2_proportion'].values
    non_responder_vals = sample_props[sample_props['stomach_pre_grouping'] == 'No-response']['C2_proportion'].values

    stat, pval = stats.mannwhitneyu(non_responder_vals, responder_vals, alternative='two-sided')
    print(f"Responders (n={len(responder_vals)}): mean={np.mean(responder_vals):.2f}%")
    print(f"Non-Responders (n={len(non_responder_vals)}): mean={np.mean(non_responder_vals):.2f}%")
    print(f"Mann-Whitney p={pval:.4f}")

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

    data = [responder_vals, non_responder_vals]
    # ONE BOX (2026-09-16, the author's fifth reading: "box plots should
    # look alike throughout"): cnsfig.boxes.draw_boxes - the 2H/I box with
    # 0.79 mm open fliers and a black median (the per-group median colours
    # of the submitted panel are the one thing this gives up). The bracket
    # keeps its vertices (line at 1.15 x y_max, arms 5% of that) and prints
    # through p_text_kw, the paper's one P-label owner: 0.0571 is "P = 0.06"
    # here as on S3 D, where it was "P = 0.057" from a local .3f (R1.3c asked
    # for exact values; the author's ruling of 2026-09-14 is two decimals at
    # or above 0.05). Its ink sits 0.4 mm above the line.
    bp = draw_boxes(ax, data, [1, 2],
                    [COLORS['Responsed'], COLORS['No-response']], width=0.6)

    y_max = max(np.max(responder_vals), np.max(non_responder_vals))
    _, p_text, _ = bracket(fig, ax, 1, 2, y_max, y_max, pval, kind="pair",
                           lift=0.15, arm=0.05 * 1.15)

    ax.set_title('Epi_CEACAM5/6')
    ax.set_ylabel('Proportion (%)')
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['Pre-R', 'Pre-NR'])
    ax.set_ylim(0, y_max * 1.40)
    ylim_above(ax, p_text)
    assert_no_points(ax, bp)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_CELL)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    style.save_panel(fig, Path(OUTPUT_DIR) / "c2_proportion_pre_boxplot")
    print(f"  Saved: c2_proportion_pre_boxplot.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
