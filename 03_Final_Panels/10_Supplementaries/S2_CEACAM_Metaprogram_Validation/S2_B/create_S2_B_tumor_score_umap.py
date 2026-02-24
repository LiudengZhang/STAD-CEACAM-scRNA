#!/usr/bin/env python3
"""
S2_B: Tumor Score UMAP on epithelial cells
4x scaling method for crisp text rendering

Tumor score (malignant_score - non_malignant_score, Seurat AddModuleScore).
Transferred via cell barcode matching (100% overlap confirmed).
"""

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

# Central config
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import EPITHELIAL_H5AD, EPITHELIAL_TUMOR_SCORED_H5AD

# Nature Cancer specifications - 4x scaling method
DPI = 300
PANEL_WIDTH_CM = 3.5 * 4
PANEL_HEIGHT_CM = 3.5 * 4
CM_TO_INCH = 1 / 2.54
SCALE = 4

# Paths
OUTPUT_DIR = Path(__file__).parent
OUTPUT_FILE = OUTPUT_DIR / "panel_S2_B.png"


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

    # Transfer tumor_score from tumor-scored file
    print("Loading tumor scores...")
    adata_ts = sc.read_h5ad(EPITHELIAL_TUMOR_SCORED_H5AD, backed='r')
    tumor_scores = adata_ts.obs['tumor_score']
    adata.obs['tumor_score'] = tumor_scores.reindex(adata.obs.index)
    n_valid = adata.obs['tumor_score'].notna().sum()
    print(f"  Transferred tumor_score for {n_valid}/{adata.n_obs} cells")

    fig_width = PANEL_WIDTH_CM * CM_TO_INCH
    fig_height = PANEL_HEIGHT_CM * CM_TO_INCH
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    sc.pl.umap(
        adata,
        color='tumor_score',
        ax=ax,
        show=False,
        frameon=False,
        title='',
        size=3 * SCALE,
        cmap='RdBu_r',
        vmin=-3,
        vmax=3,
        colorbar_loc='right',
    )

    # Rasterize scatter dots
    for coll in ax.collections:
        coll.set_rasterized(True)

    ax.set_xlabel('', fontsize=6 * SCALE)
    ax.set_ylabel('', fontsize=6 * SCALE)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Tumor Score', fontsize=9 * SCALE, fontweight='normal')

    # Remove frame
    for spine in ax.spines.values():
        spine.set_visible(False)

    cbar = ax.collections[0].colorbar
    if cbar:
        cbar.ax.tick_params(labelsize=5 * SCALE, width=1.0 * SCALE, length=4 * SCALE)
        cbar.set_label('Score', fontsize=5 * SCALE)
        cbar.outline.set_linewidth(1.0 * SCALE)

    plt.tight_layout()

    plt.savefig(OUTPUT_FILE, dpi=DPI, facecolor='white', bbox_inches='tight')
    plt.savefig(OUTPUT_FILE.with_suffix('.svg'), format='svg', dpi=DPI, facecolor='white', bbox_inches='tight')
    print(f"Saved: {OUTPUT_FILE}")
    plt.close()


if __name__ == '__main__':
    main()
