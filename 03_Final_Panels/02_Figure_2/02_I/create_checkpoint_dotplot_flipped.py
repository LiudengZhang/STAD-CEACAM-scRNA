#!/usr/bin/env python3
"""
Panel H: Flipped Checkpoint Dotplot
- Genes on X-axis (horizontal)
- Log2 FC (NR/R) on Y-axis
- Wide format for Row 2 layout
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
from paths import *

# Nature Cancer specifications
DPI = 300
SCALE = 4
PANEL_WIDTH_CM = 11.5 * SCALE   # Wide format
PANEL_HEIGHT_CM = 5.5 * SCALE   # Shorter
CM_TO_INCH = 1 / 2.54

# Paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
# DATA_PATH - now using EPITHELIAL_H5AD from central config
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
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    print("Loading epithelial data...")
    adata = sc.read_h5ad(EPITHELIAL_H5AD)
    print(f"Loaded {adata.n_obs} cells")

    data = compute_checkpoint_statistics(adata)
    # Sort by Log2FC for better visualization
    data = data.sort_values('Merged_Log2FC', ascending=False)

    genes = data['Gene'].values
    fold_changes = data['Merged_Log2FC'].values
    p_values = data['Merged_P_Value'].values
    sizes = [assign_dot_size(p) for p in p_values]

    # FLIPPED: genes on X-axis, Log2FC on Y-axis
    x_positions = np.arange(len(genes))

    fig_width = PANEL_WIDTH_CM * CM_TO_INCH
    fig_height = PANEL_HEIGHT_CM * CM_TO_INCH
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    cmap = plt.cm.RdBu_r

    # FLIPPED scatter: x=gene positions, y=fold_changes
    scatter = ax.scatter(
        x_positions, fold_changes, s=sizes, c=fold_changes, cmap=cmap, alpha=0.8,
        edgecolors='black', linewidths=0.5, vmin=-2, vmax=2
    )

    # X-axis: genes
    ax.set_xticks(x_positions)
    ax.set_xticklabels(genes, fontsize=5 * SCALE, rotation=45, ha='right', style='italic')
    ax.set_xlim(-1, len(genes))

    # Y-axis: Log2FC
    y_min, y_max = fold_changes.min(), fold_changes.max()
    y_range = y_max - y_min
    ax.set_ylim(y_min - 0.15 * y_range, y_max + 0.15 * y_range)
    ax.set_ylabel('Log2 FC (NR/R)', fontsize=7 * SCALE)

    # Horizontal reference line at 0
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.7)
    ax.grid(True, axis='y', alpha=0.3, linestyle=':', linewidth=0.3)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.tick_params(axis='both', labelsize=6 * SCALE, width=0.5, length=4)

    # Legend for significance
    legend_elements = [
        plt.scatter([], [], s=150*SCALE, c='gray', alpha=0.6, edgecolors='black', linewidths=0.5, label='p ≤ 0.10'),
        plt.scatter([], [], s=100*SCALE, c='gray', alpha=0.6, edgecolors='black', linewidths=0.5, label='p ≤ 0.15'),
        plt.scatter([], [], s=60*SCALE, c='gray', alpha=0.6, edgecolors='black', linewidths=0.5, label='p > 0.15')
    ]

    legend = ax.legend(handles=legend_elements, title='Significance', loc='upper right',
                       frameon=True, fontsize=5*SCALE, title_fontsize=6*SCALE,
                       handletextpad=0.2, borderpad=0.4, edgecolor='black', framealpha=0.9)
    legend.get_frame().set_linewidth(0.5)

    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label('Log2 FC', fontsize=6*SCALE)
    cbar.ax.tick_params(labelsize=5*SCALE)

    plt.tight_layout()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    plt.savefig(os.path.join(OUTPUT_DIR, "checkpoint_dotplot_flipped.png"), dpi=DPI, facecolor='white', bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, "checkpoint_dotplot_flipped.svg"), format='svg', facecolor='white', bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, "checkpoint_dotplot_flipped.pdf"), dpi=DPI, facecolor='white', bbox_inches='tight')
    print(f"\nSaved to {OUTPUT_DIR}")
    plt.close()


if __name__ == '__main__':
    main()
