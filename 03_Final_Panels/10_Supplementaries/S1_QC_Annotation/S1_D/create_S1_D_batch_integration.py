#!/usr/bin/env python3
"""
S1_D: Batch Integration - Before/After Harmonization UMAPs
Shows UMAP colored by sample ID before and after batch correction.
Style matches main figure UMAPs (frameon=False, alpha=0.6, rasterized).
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import FULL_DATASET_H5AD

import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# Configuration — 4× scaling
# =============================================================================
SCALE = 4
DPI = 300
CM_TO_INCH = 1 / 2.54
PANEL_WIDTH_CM = 16 * SCALE
PANEL_HEIGHT_CM = 7 * SCALE
N_CELLS = 30000

OUTPUT_DIR = Path(__file__).parent
OUTPUT_FILE = OUTPUT_DIR / "panel_S1_D.png"

# =============================================================================
# Load and subsample data
# =============================================================================
print("Loading data (backed mode)...")
adata = sc.read_h5ad(FULL_DATASET_H5AD, backed='r')
print(f"Full dataset: {adata.n_obs} cells")

np.random.seed(42)
n_plot = min(N_CELLS, adata.n_obs)
idx = np.sort(np.random.choice(adata.n_obs, n_plot, replace=False))

print(f"Subsampling to {n_plot} cells...")

sample_ids = adata.obs['Sample ID'].iloc[idx].values
X_pca_uncorrected = adata.obsm['X_pca_uncorrected'][idx]
X_umap_corrected = adata.obsm['X_umap'][idx]

print("Computing UMAP from uncorrected PCA...")
import anndata
adata_sub = anndata.AnnData(
    X=np.zeros((n_plot, 10)),
    obsm={'X_pca': X_pca_uncorrected}
)
adata_sub.obs['Sample ID'] = sample_ids

sc.pp.neighbors(adata_sub, use_rep='X_pca', n_neighbors=15)
sc.tl.umap(adata_sub)
X_umap_uncorrected = adata_sub.obsm['X_umap']

# =============================================================================
# Color mapping
# =============================================================================
samples = np.unique(sample_ids)
n_samples = len(samples)
print(f"Plotting {n_samples} samples")

cmap = plt.cm.get_cmap('Set2', 8)

sample_to_idx = {s: i for i, s in enumerate(samples)}
colors = [cmap(sample_to_idx[s] % 8) for s in sample_ids]

# =============================================================================
# Create figure — main figure UMAP style
# =============================================================================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
    'font.size': 5 * SCALE,
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

fig, axes = plt.subplots(1, 2, figsize=(
    PANEL_WIDTH_CM * CM_TO_INCH,
    PANEL_HEIGHT_CM * CM_TO_INCH
))
plt.subplots_adjust(wspace=0.3)

for ax, umap_coords, title in [
    (axes[0], X_umap_uncorrected, 'Before Harmonization'),
    (axes[1], X_umap_corrected, 'After Harmonization'),
]:
    ax.scatter(umap_coords[:, 0], umap_coords[:, 1],
               c=colors, s=8, alpha=0.6, rasterized=True)
    ax.set_title(title, fontsize=5.5 * SCALE)
    ax.set_xlabel('UMAP1', fontsize=5 * SCALE)
    ax.set_ylabel('UMAP2', fontsize=5 * SCALE)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

plt.savefig(OUTPUT_FILE, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.savefig(OUTPUT_FILE.with_suffix('.svg'), format='svg', dpi=DPI,
            bbox_inches='tight', facecolor='white')
plt.close()

print(f"Saved: {OUTPUT_FILE}")
print("Done!")
