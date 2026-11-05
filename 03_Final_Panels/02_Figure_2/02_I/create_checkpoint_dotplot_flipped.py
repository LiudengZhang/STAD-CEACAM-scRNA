#!/usr/bin/env python3
"""
Figure 2J: immune-checkpoint gene expression in pre-treatment stomach
epithelium, non-responders versus responders.

Each dot is one gene: height and colour are the log2 fold change (NR/R) of
the sample-mean expression; area is the P band of a two-sided Mann-Whitney U
test over samples (<= 0.10, <= 0.15, > 0.15).

Reads data/checkpoint_stats.csv beside this script; with --recompute, computes
it from the epithelial h5ad (EPITHELIAL_H5AD).
Writes checkpoint_dotplot_flipped.[svg|pdf|png] beside this script.

Run: python create_checkpoint_dotplot_flipped.py [--recompute]

Margins are set before the colour bar is made, because matplotlib sizes the
bar from the axes box.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import EPITHELIAL_H5AD                          # noqa: E402
from shared.expression import assert_x_stands_in_for_raw   # noqa: E402
from cnsfig import Panel, style  # noqa: E402
from cnsfig import cache                                  # noqa: E402

SCALE = 4                           # sets MARK relative to the tick size
SMALL_PT = 5.0

PANEL = Panel(2, "J", __file__, small_pt=SMALL_PT)
MARGIN = dict(left=6.5, right=0.8, top=2.5, bottom=12.0)

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


# Dot area for each P band, relative to the largest dot (2.5 mm across).
# scatter's s is the diameter squared in points.
DOT_D_MAX_MM = 2.5
DOT_S_MAX = (DOT_D_MAX_MM * style.PT_PER_MM) ** 2
DOT_AREA_RATIO = {"<=0.10": 1.0, "<=0.15": 0.55, ">0.15": 0.15}


def assign_dot_size(p_value):
    if p_value <= 0.10:
        return DOT_S_MAX * DOT_AREA_RATIO["<=0.10"]
    elif p_value <= 0.15:
        return DOT_S_MAX * DOT_AREA_RATIO["<=0.15"]
    else:
        return DOT_S_MAX * DOT_AREA_RATIO[">0.15"]


def compute_checkpoint_statistics(adata):
    pre_mask = adata.obs['Treatment phase'] == 'Pre'
    valid_mask = adata.obs['stomach_pre_grouping'].isin(['Responsed', 'No-response'])
    adata_filtered = adata[pre_mask & valid_mask].copy()

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
    style.apply(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2

    def compute():
        adata = sc.read_h5ad(EPITHELIAL_H5AD)
        # Raises if this file has no .raw and its .X is not known to equal it.
        assert_x_stands_in_for_raw(adata, EPITHELIAL_H5AD)
        return compute_checkpoint_statistics(adata)

    data = cache.table(PANEL.dir, "checkpoint_stats", compute)
    data = data.sort_values('Merged_Log2FC', ascending=False)

    genes = data['Gene'].values
    fold_changes = data['Merged_Log2FC'].values
    p_values = data['Merged_P_Value'].values
    sizes = [assign_dot_size(p) for p in p_values]

    # Genes on the x axis, Log2FC on the y
    x_positions = np.arange(len(genes))

    fig, ax = style.subplots_mm(PANEL.w_mm, PANEL.h_mm)
    # Before the colour bar (see header).
    style.margins_mm(fig, **MARGIN)

    cmap = plt.cm.RdBu      # higher in NR (positive log2 FC) blue, higher in R red

    scatter = ax.scatter(
        x_positions, fold_changes, s=sizes, c=fold_changes, cmap=cmap, alpha=0.8,
        edgecolors='black', linewidths=style.EDGE_PT, vmin=-2, vmax=2
    )

    ax.set_xticks(x_positions)
    # Steep rotation: at 45 degrees the 29 labels, 3.3 mm apart, overlap.
    ax.set_xticklabels(genes, rotation=80, ha='right', style='italic',
                       rotation_mode='anchor')
    ax.set_xlim(-1, len(genes))

    y_min, y_max = fold_changes.min(), fold_changes.max()
    y_range = y_max - y_min
    ax.set_ylim(y_min - 0.15 * y_range, y_max + 0.15 * y_range)
    ax.set_ylabel('Log2 FC (NR/R)')

    # Zero line only, no grid.
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=style.RULE_PT, alpha=0.7)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Significance key. The style sets legend.markerscale, which scales a
    # scatter key's area by its square; divide it out so each key matches its dot.
    key = 1.0 / plt.rcParams["legend.markerscale"] ** 2
    legend_elements = [
        plt.scatter([], [], s=assign_dot_size(0.10)*key, c='gray', alpha=0.6, edgecolors='black', linewidths=style.EDGE_PT, label='P ≤ 0.10'),
        plt.scatter([], [], s=assign_dot_size(0.15)*key, c='gray', alpha=0.6, edgecolors='black', linewidths=style.EDGE_PT, label='P ≤ 0.15'),
        plt.scatter([], [], s=assign_dot_size(0.50)*key, c='gray', alpha=0.6, edgecolors='black', linewidths=style.EDGE_PT, label='P > 0.15')
    ]

    # Extra spacing so the keys do not touch their labels.
    legend = ax.legend(handles=legend_elements, title='Significance', loc='upper right',
                       frameon=True, handletextpad=0.8, borderpad=0.5,
                       labelspacing=0.6, edgecolor='black', framealpha=0.9)
    legend.get_frame().set_linewidth(style.RULE_PT)

    cbar = plt.colorbar(scatter, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label('Log2 FC')

    os.makedirs(PANEL.dir, exist_ok=True)
    PANEL.finish(fig, 'checkpoint_dotplot_flipped')


if __name__ == '__main__':
    main()
