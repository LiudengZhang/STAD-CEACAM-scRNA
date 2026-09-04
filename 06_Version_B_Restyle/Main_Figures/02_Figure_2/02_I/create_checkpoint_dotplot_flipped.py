#!/usr/bin/env python3
"""
Figure 2, printed panel J, RESTYLED (Version B) - flipped checkpoint dotplot.

- Genes on X-axis (horizontal)
- Log2 FC (NR/R) on Y-axis
- Wide format for Row 2 layout

Version A is
`03_Final_Panels/02_Figure_2/02_I/create_checkpoint_dotplot_flipped.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md
allows: the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every gene, every filter, every statistic, every colour limit and every string
is Version A's. The drawing code is the same code.

  printed panel  Figure 2 J     (PROVENANCE.csv; NOT inferred from "02_I")
  printed rect   104.9 x 55.1 mm    (panel_rects.csv)
  Version B box  140.0 x 70.0 mm

MARK
    Version A drew at SCALE = 4 (11.5 x 5.5 cm x 4 = 460 x 220 mm) and its
    smallest body type is the gene tick labels and the legend, both at
    `5 * SCALE`. So

        SCALE = 4, SMALL_PT = 5
        MARK  = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 20 = 0.35
        AREA  = MARK ** 2 = 0.1225

    The three significance dot areas (`150 * SCALE`, `100 * SCALE`,
    `60 * SCALE`) take AREA; the dot edge width, the zero rule, the y grid and
    the legend frame take MARK.

    Version A's spine linewidths and `tick_params(width=..., length=...)` are
    NOT carried over: those are axes furniture, cnsplots has its own settings
    for them, and following the library rather than rescaling the old numbers
    is the standard-methods rule. `frameon=True` and the frame colour on the
    legend ARE carried over - they are Version A's explicit choice and are not
    type.

LEGEND KEYS
    cnsplots sets `legend.markerscale = 0.5`, which a scatter handler applies as
    an *area* factor of 0.25. Version A had no such reduction, so each key is
    divided by `markerscale ** 2` and prints at exactly the area of the dot it
    stands for. Same correction, and same reason, as the stage-2 exemplar S8_F.

LAYOUT ORDER
    `style.margins_mm` is called before the colorbar is made, not after. A
    colorbar built with `plt.colorbar(..., ax=ax)` steals its space out of the
    axes' *current* position; calling subplots_adjust afterwards would move the
    axes back over the colorbar. Version A got away with `tight_layout()` last
    because tight_layout knows about the colorbar's axes. Nothing drawn moves;
    only the order of two layout calls.

    The panel grew from 104.9 x 55.1 mm to 140 x 70 mm: twenty-nine
    45-degree gene labels at 7 pt need about 3.4 mm of column pitch and 9 mm of
    depth, and the largest dot is 3.0 mm across.
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
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import *                                       # noqa: E402,F403
import panel_style_cns as style                           # noqa: E402

SCALE = 4                           # Version A's canvas multiplier
SMALL_PT = 5.0                      # Version A's smallest body type

PRINTED_MM = (104.9, 55.1)          # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 140.0, 70.0
MARGIN = dict(left=11.0, right=11.0, top=3.0, bottom=12.0)

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
    family = style.apply()
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
    ax.set_xticklabels(genes, rotation=45, ha='right', style='italic')
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
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    style.save_panel(fig, Path(OUTPUT_DIR) / "checkpoint_dotplot_flipped")
    print(f"\nSaved: {Path(OUTPUT_DIR) / 'checkpoint_dotplot_flipped'}"
          f".[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
