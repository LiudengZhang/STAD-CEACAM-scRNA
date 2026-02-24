#!/usr/bin/env python3
"""
Panel D: TNF, IL1B, IL6, IL1A gene expression UMAPs for MoMac cells.
4× scaling method for Nature Cancer.
"""

import scanpy as sc
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *

BASE_DIR = Path(__file__).parent

# 4× scaling
DPI = 300
SCALE = 4
CM_TO_INCH = 1 / 2.54
PANEL_WIDTH_CM = 5.0 * SCALE
PANEL_HEIGHT_CM = 4.2 * SCALE

GENES = ['TNF', 'IL1B', 'IL6', 'IL1A']


def create_gene_umap(adata, gene, ax):
    if gene not in adata.var_names and (adata.raw is None or gene not in adata.raw.var_names):
        ax.text(0.5, 0.5, f'{gene}\n(not found)', ha='center', va='center', transform=ax.transAxes)
        ax.set_xticks([]); ax.set_yticks([])
        return

    umap = adata.obsm['X_umap']
    if adata.raw is not None and gene in adata.raw.var_names:
        gene_idx = list(adata.raw.var_names).index(gene)
        expression = adata.raw.X[:, gene_idx]
    else:
        gene_idx = list(adata.var_names).index(gene)
        expression = adata.X[:, gene_idx]

    if hasattr(expression, 'toarray'):
        expression = expression.toarray().flatten()
    else:
        expression = np.array(expression).flatten()

    scatter = ax.scatter(umap[:, 0], umap[:, 1], c=expression, cmap='Reds',
                         s=1, alpha=0.8,
                         vmin=np.percentile(expression, 2),
                         vmax=np.percentile(expression, 98))

    ax.set_title(gene, fontsize=7 * SCALE, style='italic')
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(''); ax.set_ylabel('')

    cbar = plt.colorbar(scatter, ax=ax, shrink=0.8, pad=0.02)
    cbar.ax.tick_params(labelsize=5 * SCALE, width=1.0, length=3)
    cbar.outline.set_linewidth(1.0)


def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    print("Loading MoMac data...")
    adata = sc.read_h5ad(MOMAC_H5AD)
    print(f"  Loaded {adata.n_obs} cells")

    if 'X_umap' not in adata.obsm:
        print("  Computing UMAP...")
        sc.pp.neighbors(adata, use_rep='X_pca')
        sc.tl.umap(adata)

    fig, axes = plt.subplots(2, 2, figsize=(PANEL_WIDTH_CM * CM_TO_INCH, PANEL_HEIGHT_CM * CM_TO_INCH))
    axes = axes.flatten()

    for idx, gene in enumerate(GENES):
        create_gene_umap(adata, gene, axes[idx])

    plt.tight_layout()

    output = BASE_DIR / 'tnf_il1b_il6_il1a_cytokines.png'
    plt.savefig(output, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.savefig(output.with_suffix('.svg'), bbox_inches='tight', facecolor='white')
    plt.savefig(output.with_suffix('.pdf'), format='pdf', dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {output}")


if __name__ == "__main__":
    main()
