#!/usr/bin/env python3
# REVISED FOR CIR-26-0753-ET, reviewer 1 point R1.3c:
# the test is two-sided and the annotation reports the exact P value.
# Original one-tailed version: Round_5/03_Final_Panels/02_Figure_2/02_J/create_ceacam6_boxplot.py
"""
Figure 2, printed panel K (CEACAM6 half), RESTYLED (Version B).

Version A is
`03_Final_Panels/02_Figure_2/02_J/create_ceacam6_boxplot.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md
allows: the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every value read, every filter, every statistic and every string is Version A's.
The drawing code is the same code.

  printed panel  Figure 2 K     (PROVENANCE.csv; NOT inferred from "02_J")
                 Printed K is ONE letter over TWO drawings: 02_J (CEACAM6,
                 this file) and 02_K (CEACAM5). They stay two files, as in
                 Version A; the assembler places them side by side.
  printed rect   138.8 x 20.5 mm for the pair (panel_rects.csv), so roughly
                 69 x 20.5 mm for this drawing alone
  Version B box  42.0 x 39.0 mm

MARK
    Version A drew at SCALE = 4 (`panel_figsize(3.2, 3.0)` = 128 x 120 mm) and
    its smallest body type is the 6 pt nominal on the y label, the tick labels
    and the bracket P value - `6 * SCALE`. So

        SCALE = 4, SMALL_PT = 6
        MARK  = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 24 = 0.29167
        AREA  = MARK ** 2

    Every non-type length - box, whisker, cap, median and bracket line widths,
    the flier marker and its edge - is multiplied by MARK, so each keeps the
    size it had relative to the type.

    Version A's spine linewidths and `tick_params(width=..., length=...)` are
    NOT carried over: those are axes furniture, cnsplots has its own settings
    for them (axes.linewidth 0.5, tick width 0.6, length 2), and following the
    library rather than rescaling the old numbers is the standard-methods rule.
    Nothing they control is a plotted value.
"""

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import EPITHELIAL_H5AD                                    # noqa: E402
from shared.figure_config import (RESPONSE_COLORS,                   # noqa: E402
                                  RESPONSE_MEDIAN_COLORS)
import panel_style_cns as style                                      # noqa: E402

COLORS = RESPONSE_COLORS
MEDIAN_COLORS = RESPONSE_MEDIAN_COLORS

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = BASE_DIR

PRINTED_MM = (138.8, 20.5)          # published rect of printed K (both halves)
PANEL_W_MM, PANEL_H_MM = 42.0, 39.0
MARGIN = dict(left=9.0, right=2.5, top=5.5, bottom=6.0)

SCALE = 4                           # Version A's canvas multiplier
SMALL_PT = 6.0                      # Version A's smallest body type


def main():
    family = style.apply()
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
    # inputs is the scaled matrix left by the double normalisation of
    # 2025-07-30 and carries NaN rows. The earlier version fell back to it when
    # CEACAM6 was missing from .raw, which would have drawn a panel instead of
    # stopping - that is not what the guard below does. The clean deposit
    # promotes .raw.X to .X and carries no .raw at all, so a file without one
    # holds the same log1p numbers in .X.
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

    # Create boxplot - line widths held at their Version A size relative to type
    bp = ax.boxplot(
        data,
        positions=positions,
        widths=0.6,
        patch_artist=True,
        showfliers=True,
        boxprops=dict(linewidth=0.5 * MARK),
        whiskerprops=dict(color='black', linewidth=0.5 * MARK),
        capprops=dict(color='black', linewidth=0.5 * MARK),
        flierprops=dict(marker='o', markerfacecolor='white',
                        markersize=4 * MARK,
                        linestyle='none', markeredgecolor='black',
                        markeredgewidth=0.5 * MARK)
    )

    # Style boxes
    bp['boxes'][0].set_facecolor(COLORS['Responsed'])
    bp['boxes'][0].set_edgecolor('black')
    bp['boxes'][1].set_facecolor(COLORS['No-response'])
    bp['boxes'][1].set_edgecolor('black')

    # Color median lines
    bp['medians'][0].set_color(MEDIAN_COLORS['Responsed'])
    bp['medians'][0].set_linewidth(0.8 * MARK)
    bp['medians'][1].set_color(MEDIAN_COLORS['No-response'])
    bp['medians'][1].set_linewidth(0.8 * MARK)

    # Add significance bracket
    y_max = max(np.max(responder_vals), np.max(non_responder_vals))
    y_bracket = y_max * 1.15
    ax.plot([1, 1, 2, 2], [y_bracket, y_bracket*1.05, y_bracket*1.05, y_bracket],
            'k-', linewidth=0.5 * MARK)
    pval_text = f'P = {pval:.3f}' if pval >= 0.001 else 'P < 0.001'
    ax.text(1.5, y_bracket*1.08, pval_text, ha='center', va='bottom',
            fontsize=style.tick_pt())

    # Labels - sizes now come from cnsplots' rcParams
    ax.set_title('CEACAM6')
    ax.set_ylabel('Expression')
    ax.set_xticks(positions)
    ax.set_xticklabels(['R', 'NR'])
    ax.set_ylim(0, y_max * 1.40)

    # Styling (no grid per user preference)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Margins in millimetres of paper, not fractions of a 4x canvas
    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    style.save_panel(fig, Path(OUTPUT_DIR) / "ceacam6_pre_boxplot")
    print(f"\nSaved: {Path(OUTPUT_DIR) / 'ceacam6_pre_boxplot'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
