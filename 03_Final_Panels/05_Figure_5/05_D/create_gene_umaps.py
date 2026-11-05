#!/usr/bin/env python3
"""Figure 5C: TNF, IL1B, IL6 and IL1A expression on the MoMac UMAP, 2 x 2,
colour limits at the 2nd and 98th percentiles. Expression is read through
shared.expression.expression_source and checked to be log1p(CP10k).

Input:  MOMAC_H5AD (00_Config/paths.py)
Output: tnf_il1b_il6_il1a_cytokines.{svg,pdf,png} in this directory
Run:    python create_gene_umaps.py
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
from matrix_check import assert_log1p_cp10k
from shared.expression import expression_source
from cnsfig import Panel, style
from cnsfig.layout import umap_grid_mm

SCALE = 4                       # marker-size scaling constants
SMALL_PT = 5.0

PANEL = Panel(5, "C", __file__, small_pt=SMALL_PT)
#: Map size and grid offsets in mm, shared with panel B so the maps match.
MAP_MM, GRID_LEFT_MM, GRID_TOP_MM = 12.8, 4.3, 1.0

# Set in main() once the style is applied; used by create_gene_umap().
AREA = None

GENES = ['TNF', 'IL1B', 'IL6', 'IL1A']


def create_gene_umap(adata, src, gene, ax, cax):
    # A gene missing from the matrix is an error.
    if gene not in src.var_names:
        raise RuntimeError(
            f"{gene} is declared in GENES and is not in "
            f"{Path(MOMAC_H5AD).name}'s expression matrix. Refusing to draw a "
            "placeholder into the panel.")

    umap = adata.obsm['X_umap']
    gene_idx = list(src.var_names).index(gene)
    expression = src.X[:, gene_idx]

    if hasattr(expression, 'toarray'):
        expression = expression.toarray().flatten()
    else:
        expression = np.array(expression).flatten()

    scatter = ax.scatter(umap[:, 0], umap[:, 1], c=expression, cmap='Reds',
                         s=1 * AREA, alpha=0.8,
                         vmin=np.percentile(expression, 2),
                         vmax=np.percentile(expression, 98))
    scatter.set_rasterized(True)

    ax.set_title(gene, style='italic', pad=2.0)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(''); ax.set_ylabel('')

    cbar = plt.colorbar(scatter, cax=cax)
    cbar.ax.tick_params(labelsize=style.tick_pt(), width=style.RULE_PT, length=1.5, pad=1.0)
    cbar.outline.set_linewidth(style.RULE_PT)


def main():
    global AREA
    PANEL.apply()
    AREA = (style.tick_pt() / (SMALL_PT * SCALE)) ** 2

    adata = sc.read_h5ad(MOMAC_H5AD)
    src = expression_source(adata, MOMAC_H5AD)
    assert_log1p_cp10k(src.X, f"{Path(MOMAC_H5AD).name} expression")

    if 'X_umap' not in adata.obsm:
        sc.pp.neighbors(adata, use_rep='X_pca')
        sc.tl.umap(adata)

    fig = style.figure_mm(PANEL.w_mm, PANEL.h_mm)
    cells = umap_grid_mm(fig, n_rows=2, n_cols=2, map_mm=MAP_MM,
                         left_mm=GRID_LEFT_MM, top_mm=GRID_TOP_MM)
    axes = [ax for ax, _ in cells]
    caxes = [cax for _, cax in cells]

    for idx, gene in enumerate(GENES):
        create_gene_umap(adata, src, gene, axes[idx], caxes[idx])

    PANEL.finish(fig, 'tnf_il1b_il6_il1a_cytokines')


if __name__ == "__main__":
    main()
