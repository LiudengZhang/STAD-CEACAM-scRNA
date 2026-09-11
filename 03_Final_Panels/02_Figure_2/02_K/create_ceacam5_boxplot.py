#!/usr/bin/env python3
# REVISED FOR CIR-26-0753-ET, reviewer 1 point R1.3c:
# the test is two-sided and the annotation reports the exact P value.
# The one-tailed version this replaces is the one used in the preprint,
# https://www.biorxiv.org/content/10.64898/2026.03.05.708917
"""
Figure 2 panel K, CEACAM5 box - sample-mean CEACAM5 expression in
pre-treatment stomach epithelium, responders versus non-responders.

Printed K is one letter over two drawings, 02_J (CEACAM6) and 02_K (CEACAM5, this file). Each is
drawn at its own printed sub-box, read from 03_Final_Panels/slot_subrects.csv
through 00_Config/slots.py, so the type size set here is the type size printed.
Margins are measured from the rendered ink rather than typed.

The printed letter K sits to the left of both sub-boxes, so neither of them
reserves a corner for it.

  printed panel  Figure 2 K, box 2   (PROVENANCE.csv; NOT inferred from "02_K")

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

PANEL_LETTER = "K"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(2, PANEL_LETTER, sub=2)

SCALE = 4                           # the earlier canvas multiplier
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

    # The working inputs keep the log1p matrix in .raw; the clean deposit
    # promotes it to .X and carries no .raw, so a file without .raw already
    # holds the same numbers in .X.
    src = adata_filtered.raw if adata_filtered.raw is not None else adata_filtered
    idx = src.var_names.get_loc('CEACAM5')
    expr = src.X[:, idx].toarray().flatten()
    adata_filtered.obs['CEACAM5'] = expr

    sample_means = adata_filtered.obs.groupby('sample', observed=True).agg({
        'CEACAM5': 'mean',
        'stomach_pre_grouping': 'first'
    }).reset_index()

    responder_vals = sample_means[sample_means['stomach_pre_grouping'] == 'Responsed']['CEACAM5'].values
    non_responder_vals = sample_means[sample_means['stomach_pre_grouping'] == 'No-response']['CEACAM5'].values

    stat, pval = stats.mannwhitneyu(non_responder_vals, responder_vals, alternative='two-sided')
    print(f"Responders (n={len(responder_vals)}): mean={np.mean(responder_vals):.3f}")
    print(f"Non-Responders (n={len(non_responder_vals)}): mean={np.mean(non_responder_vals):.3f}")
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
    pval_text = f'P = {pval:.3f}' if pval >= 0.001 else 'P < 0.001'
    ax.text(1.5, y_bracket*1.08, pval_text, ha='center', va='bottom',
            fontsize=style.tick_pt())

    # The gene symbol is set in italic, as the shipped page sets it and as
    # the rest of this figure sets one. Style only; the string is unchanged.
    ax.set_title('CEACAM5', fontstyle='italic')
    ax.set_ylabel('Expression')
    ax.set_xticks([1, 2])
    # The two response labels stand 5.5 mm apart on this box and set 5.52 and
    # 7.05 mm at 6 pt, so drawn horizontally they run into one another - by
    # 1.08 mm on the narrowest of the four boxes of panels K and L. They are
    # set at 45 degrees instead, which is what the printed page's own narrow
    # panels do: the same two strings, the same tick positions, a different
    # angle.
    ax.set_xticklabels(['Pre-R', 'Pre-NR'],
                       rotation=45, ha='right')
    ax.set_ylim(0, y_max * 1.40)
    # Three gradations rather than five: on a 19.3 mm box the fourth
    # prints against the bracket's P value. The scale, the limits and
    # the data are untouched; only how finely the ruler is marked
    # changes.
    ax.locator_params(axis='y', nbins=3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    style.fit_margins(fig, pad_mm=0.6, reserve_letter=False)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    style.save_panel(fig, Path(OUTPUT_DIR) / "ceacam5_pre_boxplot")
    print(f"Saved: {Path(OUTPUT_DIR) / 'ceacam5_pre_boxplot'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
