#!/usr/bin/env python3
"""
Panel 2D: CEACAM6 expression UMAP on epithelial cells
4× scaling method for crisp text rendering
"""

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import os

# Central config
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *
from shared.figure_config import use_panel_style

# Nature Cancer specifications - 4× scaling method
DPI = 300
PANEL_WIDTH_CM = 3.5 * 4   # 14 cm electronic → 3.5 cm print
PANEL_HEIGHT_CM = 3.5 * 4  # 14 cm electronic → 3.5 cm print
CM_TO_INCH = 1 / 2.54
SCALE = 4

# Paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
# DATA_PATH - now using EPITHELIAL_H5AD from central config
OUTPUT_DIR = BASE_DIR

def main():
    use_panel_style(font_pt=7)

    print("Loading data...")
    adata = sc.read_h5ad(EPITHELIAL_H5AD)
    print(f"Loaded {adata.n_obs} cells")

    gene = 'CEACAM6'
    # The working inputs keep the log1p matrix in .raw; the clean deposit
    # promotes it to .X and carries no .raw, so a file without .raw already
    # holds the same numbers in .X.
    src = adata.raw if adata.raw is not None else adata
    idx = src.var_names.get_loc(gene)
    expr = src.X[:, idx]
    if hasattr(expr, 'toarray'):
        expr = expr.toarray().flatten()
    else:
        expr = np.array(expr).flatten()
    adata.obs['CEACAM6_expr'] = expr

    fig_width = PANEL_WIDTH_CM * CM_TO_INCH
    fig_height = PANEL_HEIGHT_CM * CM_TO_INCH
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    sc.pl.umap(
        adata,
        color='CEACAM6_expr',
        ax=ax,
        show=False,
        frameon=False,
        title='',
        size=3 * SCALE,
        cmap='Reds',
        vmin=2,
        vmax=6,
        colorbar_loc='right',
    )

    # Rasterize scatter dots (keeps axes/legend as vectors, dots as embedded raster)
    for coll in ax.collections:
        coll.set_rasterized(True)

    ax.set_xlabel('', fontsize=6 * SCALE)
    ax.set_ylabel('', fontsize=6 * SCALE)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('CEACAM6', fontsize=9 * SCALE, fontweight='normal')

    # Remove frame
    for spine in ax.spines.values():
        spine.set_visible(False)

    cbar = ax.collections[0].colorbar
    if cbar:
        cbar.ax.tick_params(labelsize=5 * SCALE, width=1.0 * SCALE, length=4 * SCALE)
        cbar.set_label('Expression', fontsize=5 * SCALE)
        cbar.outline.set_linewidth(1.0 * SCALE)

    plt.tight_layout()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    plt.savefig(os.path.join(OUTPUT_DIR, "ceacam6_expression_umap.png"), dpi=DPI, facecolor='white', bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, "ceacam6_expression_umap.svg"), format='svg', dpi=DPI, facecolor='white', bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, "ceacam6_expression_umap.pdf"), dpi=DPI, facecolor='white', bbox_inches='tight')
    print(f"Saved to {OUTPUT_DIR}")
    plt.close()

if __name__ == '__main__':
    main()
