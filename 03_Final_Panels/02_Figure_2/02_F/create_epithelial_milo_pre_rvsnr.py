#!/usr/bin/env python3
"""
Panel 2J: Milo differential abundance analysis for epithelial cells
4× scaling method for crisp text rendering
"""

import scanpy as sc
import pertpy as pt
import mudata as mu
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')

# Central config
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *

# Nature Cancer specifications - 4× scaling method
DPI = 300
PANEL_WIDTH_CM = 4.5 * 4   # 18 cm electronic → 4.5 cm print
PANEL_HEIGHT_CM = 4.5 * 4  # 18 cm electronic → 4.5 cm print
CM_TO_INCH = 1 / 2.54
SCALE = 4

# Milo parameters
K_NEIGHBORS = 15
NHOOD_PROP = 0.1
SOLVER = "edger"
ALPHA = 0.1
MAX_CELLS = 8000

# Paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
# DATA_PATH - now using EPITHELIAL_H5AD from central config
OUTPUT_DIR = BASE_DIR

GROUP_COL = "stomach_pre_grouping"
GROUP_A = "Responsed"
GROUP_B = "No-response"
SAMPLE_COL = "sample"


def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    print("=" * 60)
    print("Milo Differential Abundance Analysis")
    print("=" * 60)

    print("\n[1/7] Loading epithelial data...")
    adata = sc.read_h5ad(EPITHELIAL_H5AD)
    print(f"  Loaded {adata.n_obs} cells")

    print("\n[2/7] Filtering for Pre-treatment cells...")
    pre_mask = adata.obs['Treatment phase'] == 'Pre'
    valid_mask = adata.obs[GROUP_COL].isin([GROUP_A, GROUP_B])
    adata = adata[pre_mask & valid_mask].copy()
    print(f"  Pre-treatment cells: {adata.n_obs}")

    if adata.n_obs > MAX_CELLS:
        print(f"\n[3/7] Subsampling to {MAX_CELLS} cells...")
        sc.pp.subsample(adata, n_obs=MAX_CELLS, random_state=42)
    else:
        print("\n[3/7] No subsampling needed")

    if 'X_pca' not in adata.obsm:
        print("\n[4/7] Computing PCA...")
        sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
    else:
        print("\n[4/7] Using existing PCA")

    print(f"\n[5/7] Computing neighbors (k={K_NEIGHBORS})...")
    sc.pp.neighbors(adata, n_neighbors=K_NEIGHBORS, n_pcs=50, key_added='neighbors')

    if 'X_umap' not in adata.obsm:
        print("  Computing UMAP...")
        sc.tl.umap(adata)

    print(f"\n[6/7] Running Milo analysis...")
    milo = pt.tl.Milo()
    mdata = mu.MuData({"rna": adata})
    milo.make_nhoods(mdata, neighbors_key="neighbors", prop=NHOOD_PROP)
    milo.count_nhoods(mdata, sample_col=SAMPLE_COL)
    milo.da_nhoods(mdata, design=f"~ {GROUP_COL}", solver=SOLVER)
    milo.build_nhood_graph(mdata)

    print("\n[7/7] Creating visualization...")
    fig_width = PANEL_WIDTH_CM * CM_TO_INCH
    fig_height = PANEL_HEIGHT_CM * CM_TO_INCH
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    milo.plot_nhood_graph(mdata, alpha=1.0, min_logFC=0, ax=ax, title="")

    ax.set_title('Pre: R vs NR', fontsize=9 * SCALE, fontweight='normal')
    ax.set_xlabel('UMAP1', fontsize=6 * SCALE)
    ax.set_ylabel('UMAP2', fontsize=6 * SCALE)
    ax.tick_params(axis='both', labelsize=5 * SCALE, width=1.0 * SCALE, length=4 * SCALE)

    for spine in ax.spines.values():
        spine.set_linewidth(1.0 * SCALE)

    plt.tight_layout()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    plt.savefig(os.path.join(OUTPUT_DIR, "epithelial_milo_pre_rvsnr.png"), dpi=DPI, facecolor='white', bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, "epithelial_milo_pre_rvsnr.svg"), format='svg', facecolor='white', bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, "epithelial_milo_pre_rvsnr.pdf"), dpi=DPI, facecolor='white', bbox_inches='tight')
    print(f"\n  Saved to {OUTPUT_DIR}")
    plt.close()

    print("\nDone!")


if __name__ == '__main__':
    main()
