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

    ax.set_title('Epi_CEACAM5/6')
    ax.set_ylabel('Proportion (%)')
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['Pre-R', 'Pre-NR'])
    ax.set_ylim(0, y_max * 1.40)
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
