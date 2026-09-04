#!/usr/bin/env python3
# Paths below refer to the upstream Round_4 processing pipeline, which is
# not part of this release. This script is included as a record of how the
# input was produced; it is not called by _run_all_panels.sh.
"""05_umap_activity.py - UMAP colored by meta-program activity"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent.parent
SCORES_FILE = Path("/path/to/Project_4_05232025/Round_4/01_Round_4.2_Standardized_Pipeline/01.1_meta_program_epithelial/04_Activity_Scores/epithelial_with_mp_scores.h5ad")
OUTPUT_DIR = BASE_DIR / "03_Output" / "05_Figures"
UMAP_SOURCE = Path("/path/to/Project_4_05232025/Round_4/04_Final_Panels/00_Set_Ups/00_Data/01_Major_Cell_Types/Epithelial_functionally_annotated.h5ad")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(SCORES_FILE)
mp_cols = [col for col in adata.obs.columns if col.endswith('_score') and col.startswith('MP')]

# Load UMAP from source file
adata_umap = sc.read_h5ad(UMAP_SOURCE)
adata.obsm['X_umap'] = adata_umap[adata.obs_names].obsm['X_umap']
del adata_umap

cmap = LinearSegmentedColormap.from_list('activity', ['#f0f0f0', '#fee0d2', '#fc9272', '#de2d26'])

n_mps = len(mp_cols)
n_cols = min(4, n_mps)
n_rows = (n_mps + n_cols - 1) // n_cols

fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))
axes = np.atleast_2d(axes)

umap_coords = adata.obsm['X_umap']

for idx, mp_col in enumerate(sorted(mp_cols)):
    row, col = idx // n_cols, idx % n_cols
    ax = axes[row, col] if n_rows > 1 else axes[0, col]
    scores = adata.obs[mp_col].values
    vmin, vmax = np.percentile(scores, [2, 98])
    scatter = ax.scatter(umap_coords[:, 0], umap_coords[:, 1], c=scores, s=1, cmap=cmap, vmin=vmin, vmax=vmax, alpha=0.8, rasterized=True)
    ax.set_title(mp_col.replace('_score', ''), fontsize=11, fontweight='bold')
    ax.set_xticks([])
    ax.set_yticks([])
    plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)

for idx in range(n_mps, n_rows * n_cols):
    row, col = idx // n_cols, idx % n_cols
    (axes[row, col] if n_rows > 1 else axes[0, col]).axis('off')

plt.suptitle('Meta-Program Activity on UMAP', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "05_umap_metaprogram_grid.png", dpi=200, bbox_inches='tight', facecolor='white')
plt.close()

print("Figure 5 complete")
