#!/usr/bin/env python3
"""
Panel 2C: CEACAM5 expression UMAP on epithelial cells
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
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    print("Loading data...")
    adata = sc.read_h5ad(EPITHELIAL_H5AD)
    print(f"Loaded {adata.n_obs} cells")

    gene = 'CEACAM5'
    idx = adata.raw.var_names.get_loc(gene)
    expr = adata.raw.X[:, idx]
    if hasattr(expr, 'toarray'):
        expr = expr.toarray().flatten()
    else:
        expr = np.array(expr).flatten()
    adata.obs['CEACAM5_expr'] = expr

    fig_width = PANEL_WIDTH_CM * CM_TO_INCH
    fig_height = PANEL_HEIGHT_CM * CM_TO_INCH
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    sc.pl.umap(
        adata,
        color='CEACAM5_expr',
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
    ax.set_title('CEACAM5', fontsize=9 * SCALE, fontweight='normal')

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
    plt.savefig(os.path.join(OUTPUT_DIR, "ceacam5_expression_umap.png"), dpi=DPI, facecolor='white', bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, "ceacam5_expression_umap.svg"), format='svg', dpi=DPI, facecolor='white', bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, "ceacam5_expression_umap.pdf"), dpi=DPI, facecolor='white', bbox_inches='tight')
    print(f"Saved to {OUTPUT_DIR}")
    plt.close()

if __name__ == '__main__':
    main()
