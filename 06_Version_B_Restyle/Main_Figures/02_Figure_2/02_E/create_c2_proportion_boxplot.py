#!/usr/bin/env python3
"""
Figure 2, printed panel E, RESTYLED (Version B) - C2_Epi_CEACAM5/6 cluster
proportion, pre-treatment responders versus non-responders.

Version A is
`03_Revised_Panels/Main_Figures/02_Figure_2/02_E/create_c2_proportion_boxplot.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md
allows: the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every value read, every filter, every statistic and every string is Version A's.
The drawing code is the same code.

  printed panel  Figure 2 E     (PROVENANCE.csv; the directory letter happens
                                 to agree here - it was still looked up)
  printed rect   30.9 x 26.1 mm    (panel_rects.csv)
  Version B box  42.0 x 42.0 mm

MARK
    Version A drew at SCALE = 4 (3.2 x 3.0 cm x 4 = 128 x 120 mm) and its
    smallest body type is the 6 pt nominal on the y label, the tick labels and
    the bracket P value - `6 * SCALE`. So

        SCALE = 4, SMALL_PT = 6
        MARK  = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 24 = 0.29167

    Version A's spine linewidths and `tick_params(width=..., length=...)` are
    NOT carried over: those are axes furniture, cnsplots has its own settings
    for them, and following the library rather than rescaling the old numbers
    is the standard-methods rule. Nothing they control is a plotted value.

    The panel grew from 30.9 x 26.1 mm to 42 x 42 mm: the two-line title
    ("CEACAM5/6 / Epithelial") is 8 pt now instead of printing at about 4 pt,
    and the y label plus its tick labels take 9 mm of the width before any data
    is drawn.
"""

import scanpy as sc
import numpy as np
from scipy import stats
import os

# Central config
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import *                                       # noqa: E402,F403
import panel_style_cns as style                           # noqa: E402

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

PRINTED_MM = (30.9, 26.1)           # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 42.0, 42.0
MARGIN = dict(left=10.0, right=2.5, top=8.5, bottom=6.0)

SCALE = 4                           # Version A's canvas multiplier
SMALL_PT = 6.0                      # Version A's smallest body type


def main():
    family = style.apply()
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
    print(f"\nResponders (n={len(responder_vals)}): mean={np.mean(responder_vals):.2f}%")
    print(f"Non-Responders (n={len(non_responder_vals)}): mean={np.mean(non_responder_vals):.2f}%")
    print(f"Mann-Whitney p={pval:.4f}")

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

    data = [responder_vals, non_responder_vals]
    bp = ax.boxplot(data, positions=[1, 2], widths=0.6, patch_artist=True,
                    boxprops=dict(linewidth=0.5 * MARK),
                    whiskerprops=dict(color='black', linewidth=0.5 * MARK),
                    capprops=dict(color='black', linewidth=0.5 * MARK),
                    flierprops=dict(marker='o', markerfacecolor='white',
                                    markersize=4 * MARK,
                                    markeredgecolor='black',
                                    markeredgewidth=0.5 * MARK))

    bp['boxes'][0].set_facecolor(COLORS['Responsed'])
    bp['boxes'][1].set_facecolor(COLORS['No-response'])
    bp['medians'][0].set_color(MEDIAN_COLORS['Responsed'])
    bp['medians'][0].set_linewidth(0.8 * MARK)
    bp['medians'][1].set_color(MEDIAN_COLORS['No-response'])
    bp['medians'][1].set_linewidth(0.8 * MARK)

    y_max = max(np.max(responder_vals), np.max(non_responder_vals))
    y_bracket = y_max * 1.15
    ax.plot([1, 1, 2, 2], [y_bracket, y_bracket*1.05, y_bracket*1.05, y_bracket],
            'k-', linewidth=0.5 * MARK)
    # Exact P rather than a threshold label; the test above is already
    # two-sided, so only the annotation had to change (R1.3c).
    pval_text = f'P = {pval:.3f}' if pval >= 0.001 else 'P < 0.001'
    ax.text(1.5, y_bracket*1.08, pval_text, ha='center', va='bottom',
            fontsize=style.tick_pt())

    ax.set_title('CEACAM5/6\nEpithelial')
    ax.set_ylabel('Proportion (%)')
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['R', 'NR'])
    ax.set_ylim(0, y_max * 1.40)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    style.save_panel(fig, Path(OUTPUT_DIR) / "c2_proportion_pre_boxplot")
    print(f"\nSaved: {Path(OUTPUT_DIR) / 'c2_proportion_pre_boxplot'}"
          f".[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
