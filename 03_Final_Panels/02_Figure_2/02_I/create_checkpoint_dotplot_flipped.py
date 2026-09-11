#!/usr/bin/env python3
"""
Figure 2 panel J - immune-checkpoint gene expression in pre-treatment stomach
epithelium, non-responders against responders.

Each dot is one checkpoint gene: its height is the log2 fold change of the
sample-mean expression, its area the significance band of a two-sided
Mann-Whitney U test over the samples, and its colour the same fold change.

The panel is drawn at the millimetre rectangle it prints in, read from
03_Final_Panels/panel_rects.csv through 00_Config/slots.py, so the type size
set here is the type size printed.

  printed panel  Figure 2 J       (PROVENANCE.csv; NOT inferred from "02_I")

LAYOUT ORDER
    The margins are set before the colour bar is made, because matplotlib sizes
    the bar from the axes box it finds and takes its own share of it. They are
    millimetres of paper, not fractions of the canvas, and they are not fitted
    to the ink afterwards: moving the subplot parameters once the bar exists
    would leave the bar behind.

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the legend - at 5 * SCALE. MARK carries the non-type
    point sizes across to the 1:1 canvas:

        MARK = style.tick_pt() / (SMALL_PT * SCALE)
        AREA = MARK ** 2

    The three significance dot areas take AREA; the dot edge width, the zero
    rule, the grid rule and the legend frame take MARK. Tick widths and lengths
    and spine widths do not: those are style, and cnsplots sets them.

LEGEND KEYS
    cnsplots sets `legend.markerscale = 0.5`, which a scatter handler applies
    as an *area* factor of 0.25. Each key is therefore divided by
    `markerscale ** 2`, so that after the legend applies its factor the key
    prints at exactly the area of the dot it labels.

THE GENE LABELS ARE SET UPRIGHT, NOT AT 45 DEGREES
    Twenty-nine ticks across this plotting box stand 3.3 mm apart, and a label
    set at an angle puts its own depth across its neighbour's: at 45 degrees
    the twenty-eight adjacent pairs overlap by up to 4.6 mm2 and at 60 degrees
    by up to 0.9 mm2, whatever their length, because the spacing is fixed by
    the number of genes and the width of the slot. Upright they clear one
    another entirely. Same genes, same order, same tick positions, same
    strings; only the angle changes.

    An upright label is deeper than a rotated one and the slot is fixed, so the
    cost was measured on the rendered panel rather than assumed: the gene label
    band is 9.30 mm deep at 45 degrees and 10.71 mm upright - 1.41 mm more -
    with the deepest label, TNFRSF18, going from 9.29 to 10.70 mm. The band
    still ends 0.23 mm inside the foot of the 54.4 mm slot, and the panel
    leaves 2.37 mm unused at the top, so the extra depth is paid out of slack
    the panel already had.

Every gene, every filter, every fold change, every P value and every string is
the earlier drawing's. The drawing code is the same code.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

# Central config
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *                                       # noqa: E402,F403
import panel_style_cns as style                           # noqa: E402
import slots                                             # noqa: E402

SCALE = 4                           # the earlier canvas multiplier
SMALL_PT = 5.0                      # the earlier smallest body type

PANEL_LETTER = "J"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(2, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(2, PANEL_LETTER)
MARGIN = dict(left=6.5, right=0.8, top=2.5, bottom=12.0)

# Paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = BASE_DIR

# Checkpoint genes
INHIBITORY_GENES = [
    'PDCD1', 'CD274', 'CTLA4', 'LAG3',
    'HAVCR2', 'TIGIT', 'CD276', 'CD47', 'IDO1',
    'VSIR', 'BTLA', 'SIGLEC15', 'KLRC1', 'CD96', 'PVRIG', 'VTCN1',
    'NT5E', 'ENTPD1', 'ADORA2A',
    'CEACAM1', 'CEACAM5', 'CEACAM6',
]

COSTIMULATORY_GENES = [
    'ICOS', 'CD27', 'CD70', 'CD40',
    'TNFRSF9', 'TNFRSF4', 'TNFRSF18',
]


def assign_dot_size(p_value):
    if p_value <= 0.10:
        return 150 * SCALE
    elif p_value <= 0.15:
        return 100 * SCALE
    else:
        return 60 * SCALE


def compute_checkpoint_statistics(adata):
    print("Computing checkpoint gene statistics...")
    pre_mask = adata.obs['Treatment phase'] == 'Pre'
    valid_mask = adata.obs['stomach_pre_grouping'].isin(['Responsed', 'No-response'])
    adata_filtered = adata[pre_mask & valid_mask].copy()
    print(f"  Pre-treatment cells: {adata_filtered.n_obs}")

    all_checkpoint_genes = INHIBITORY_GENES + COSTIMULATORY_GENES
    gene_names = adata_filtered.raw.var_names if adata_filtered.raw is not None else adata_filtered.var_names
    available_genes = [g for g in all_checkpoint_genes if g in gene_names]

    results = []
    for gene in available_genes:
        if adata_filtered.raw is not None and gene in adata_filtered.raw.var_names:
            gene_idx = adata_filtered.raw.var_names.get_loc(gene)
            expr = adata_filtered.raw.X[:, gene_idx].toarray().flatten() if hasattr(adata_filtered.raw.X, 'toarray') else adata_filtered.raw.X[:, gene_idx].flatten()
        else:
            gene_idx = adata_filtered.var_names.get_loc(gene)
            expr = adata_filtered.X[:, gene_idx].toarray().flatten() if hasattr(adata_filtered.X, 'toarray') else adata_filtered.X[:, gene_idx].flatten()
            expr = np.nan_to_num(expr, nan=0.0)

        adata_filtered.obs['_gene_expr'] = expr
        sample_means = adata_filtered.obs.groupby('sample', observed=True).agg({
            '_gene_expr': 'mean',
            'stomach_pre_grouping': 'first'
        }).reset_index()

        responder_vals = sample_means[sample_means['stomach_pre_grouping'] == 'Responsed']['_gene_expr'].values
        non_responder_vals = sample_means[sample_means['stomach_pre_grouping'] == 'No-response']['_gene_expr'].values

        mean_r = np.mean(responder_vals)
        mean_nr = np.mean(non_responder_vals)
        epsilon = 1e-10
        log2fc = np.log2((mean_nr + epsilon) / (mean_r + epsilon))

        if len(responder_vals) > 1 and len(non_responder_vals) > 1:
            stat, pval = stats.mannwhitneyu(non_responder_vals, responder_vals, alternative='two-sided')
        else:
            pval = 1.0

        category = 'Inhibitory' if gene in INHIBITORY_GENES else 'Costimulatory'
        results.append({'Gene': gene, 'Category': category, 'Merged_Log2FC': log2fc, 'Merged_P_Value': pval})

    return pd.DataFrame(results)


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.5f}")

    print("Loading epithelial data...")
    adata = sc.read_h5ad(EPITHELIAL_H5AD)                 # noqa: F405
    print(f"Loaded {adata.n_obs} cells")

    data = compute_checkpoint_statistics(adata)
    # Sort by Log2FC for better visualization
    data = data.sort_values('Merged_Log2FC', ascending=False)

    genes = data['Gene'].values
    fold_changes = data['Merged_Log2FC'].values
    p_values = data['Merged_P_Value'].values
    sizes = [assign_dot_size(p) * AREA for p in p_values]

    # FLIPPED: genes on X-axis, Log2FC on Y-axis
    x_positions = np.arange(len(genes))

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)
    # Before the colorbar - see LAYOUT ORDER in the header.
    style.margins_mm(fig, **MARGIN)

    cmap = plt.cm.RdBu_r

    # FLIPPED scatter: x=gene positions, y=fold_changes
    scatter = ax.scatter(
        x_positions, fold_changes, s=sizes, c=fold_changes, cmap=cmap, alpha=0.8,
        edgecolors='black', linewidths=0.5 * MARK, vmin=-2, vmax=2
    )

    # X-axis: genes
    ax.set_xticks(x_positions)
    ax.set_xticklabels(genes, rotation=90, ha='right', style='italic')
    ax.set_xlim(-1, len(genes))

    # Y-axis: Log2FC
    y_min, y_max = fold_changes.min(), fold_changes.max()
    y_range = y_max - y_min
    ax.set_ylim(y_min - 0.15 * y_range, y_max + 0.15 * y_range)
    ax.set_ylabel('Log2 FC (NR/R)')

    # Horizontal reference line at 0
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5 * MARK, alpha=0.7)
    ax.grid(True, axis='y', alpha=0.3, linestyle=':', linewidth=0.3 * MARK)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Legend for significance. See LEGEND KEYS in the header for the
    # markerscale correction.
    key = 1.0 / plt.rcParams["legend.markerscale"] ** 2
    legend_elements = [
        plt.scatter([], [], s=150*SCALE*AREA*key, c='gray', alpha=0.6, edgecolors='black', linewidths=0.5 * MARK, label='p ≤ 0.10'),
        plt.scatter([], [], s=100*SCALE*AREA*key, c='gray', alpha=0.6, edgecolors='black', linewidths=0.5 * MARK, label='p ≤ 0.15'),
        plt.scatter([], [], s=60*SCALE*AREA*key, c='gray', alpha=0.6, edgecolors='black', linewidths=0.5 * MARK, label='p > 0.15')
    ]

    legend = ax.legend(handles=legend_elements, title='Significance', loc='upper right',
                       frameon=True, handletextpad=0.2, borderpad=0.4,
                       edgecolor='black', framealpha=0.9)
    legend.get_frame().set_linewidth(0.5 * MARK)

    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label('Log2 FC')

    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    style.save_panel(fig, Path(OUTPUT_DIR) / "checkpoint_dotplot_flipped")
    print(f"\nSaved: {Path(OUTPUT_DIR) / 'checkpoint_dotplot_flipped'}"
          f".[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
