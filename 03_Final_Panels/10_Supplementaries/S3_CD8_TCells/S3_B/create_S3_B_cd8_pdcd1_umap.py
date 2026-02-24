#!/usr/bin/env python3
"""
S3_B: CD8+ T cell UMAP colored by PDCD1 (PD-1) expression.
Stomach samples only.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import TCD8_H5AD

import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib

SCALE = 4
DPI = 300
PANEL_WIDTH_CM = 8 * SCALE
PANEL_HEIGHT_CM = 7 * SCALE
OUTPUT_DIR = Path(__file__).parent
OUTPUT_PNG = OUTPUT_DIR / "panel_S3_B.png"

print("Loading CD8 data...")
adata = sc.read_h5ad(TCD8_H5AD)
adata = adata[adata.obs['Sample site'] == 'Stomach'].copy()
print(f"Stomach CD8: {adata.n_obs} cells")

plt.rcParams.update({
    'font.size': 6 * SCALE,
    'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

fig, ax = plt.subplots(figsize=(PANEL_WIDTH_CM / 2.54, PANEL_HEIGHT_CM / 2.54))

n_cells = adata.n_obs
point_size = max(0.5, min(4.0, 50000 / n_cells)) * SCALE

sc.pl.umap(
    adata,
    color='PDCD1',
    ax=ax,
    show=False,
    frameon=True,
    title='PDCD1 (PD-1) expression',
    size=point_size,
    cmap='Reds',
    vmin=0,
)

ax.set_xlabel('UMAP1', fontsize=6 * SCALE)
ax.set_ylabel('UMAP2', fontsize=6 * SCALE)
ax.tick_params(labelsize=5 * SCALE, width=0.5 * SCALE, length=3 * SCALE)
for spine in ax.spines.values():
    spine.set_linewidth(0.5 * SCALE)

plt.tight_layout()
fig.savefig(OUTPUT_PNG, dpi=DPI, bbox_inches='tight', facecolor='white')
fig.savefig(OUTPUT_PNG.with_suffix('.svg'), format='svg', bbox_inches='tight', facecolor='white')
plt.close('all')
print(f"Saved: {OUTPUT_PNG}")
