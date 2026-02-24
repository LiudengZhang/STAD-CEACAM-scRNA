#!/usr/bin/env python3
"""
S1_A: QC Metric UMAPs (4 panels)
Shows UMAP colored by n_genes_by_counts, total_counts, pct_counts_mt, doublet_score.
Subsampled to 30k cells for speed. Style matches main figure UMAPs.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import FULL_DATASET_H5AD

import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

# =============================================================================
# Configuration — 4x scaling
# =============================================================================
SCALE = 4
DPI = 300
CM_TO_INCH = 1 / 2.54
PANEL_WIDTH_CM = 24 * SCALE      # full width
PANEL_HEIGHT_CM = 5.5 * SCALE    # compact single row
N_CELLS = 30000

OUTPUT_DIR = Path(__file__).parent
OUTPUT_FILE = OUTPUT_DIR / "panel_S1_A.png"

# =============================================================================
# Load and subsample
# =============================================================================
print("Loading data (backed mode)...")
adata = sc.read_h5ad(FULL_DATASET_H5AD, backed='r')
print(f"Full dataset: {adata.n_obs} cells")

np.random.seed(42)
n_plot = min(N_CELLS, adata.n_obs)
idx = np.sort(np.random.choice(adata.n_obs, n_plot, replace=False))

print(f"Subsampling to {n_plot} cells...")
umap_coords = adata.obsm['X_umap'][idx]
n_genes = adata.obs['n_genes_by_counts'].iloc[idx].values.astype(float)
total_counts = adata.obs['total_counts'].iloc[idx].values.astype(float)
pct_mt = adata.obs['pct_counts_mt'].iloc[idx].values.astype(float)
doublet_score = adata.obs['doublet_score'].iloc[idx].values.astype(float)

# =============================================================================
# Create figure — 1x4 row
# =============================================================================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
    'font.size': 5 * SCALE,
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

fig, axes = plt.subplots(1, 4, figsize=(
    PANEL_WIDTH_CM * CM_TO_INCH,
    PANEL_HEIGHT_CM * CM_TO_INCH,
))
plt.subplots_adjust(wspace=0.35)

panels = [
    (n_genes, 'viridis', 'Genes per Cell', None, None),
    (total_counts, 'viridis', 'Total Counts', None, None),
    (pct_mt, 'viridis', 'Mitochondrial %', 0, 25),
    (doublet_score, 'viridis', 'Doublet Score', 0, None),
]

for ax, (values, cmap, title, vmin, vmax) in zip(axes, panels):
    # Shuffle plot order so high values render on top
    order = np.argsort(values)
    sc_plot = ax.scatter(
        umap_coords[order, 0], umap_coords[order, 1],
        c=values[order], cmap=cmap, s=8, alpha=0.6, rasterized=True,
        vmin=vmin, vmax=vmax,
    )
    ax.set_title(title, fontsize=5.5 * SCALE, pad=4)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Colorbar
    cbar = plt.colorbar(sc_plot, ax=ax, shrink=0.7, pad=0.02, aspect=20)
    cbar.ax.tick_params(labelsize=4 * SCALE)

# UMAP axis labels on first panel only
axes[0].set_xlabel('UMAP1', fontsize=4.5 * SCALE)
axes[0].set_ylabel('UMAP2', fontsize=4.5 * SCALE)

plt.savefig(OUTPUT_FILE, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.savefig(OUTPUT_FILE.with_suffix('.svg'), format='svg', dpi=DPI,
            bbox_inches='tight', facecolor='white')
plt.close()

print(f"Saved: {OUTPUT_FILE}")
print("Done!")
