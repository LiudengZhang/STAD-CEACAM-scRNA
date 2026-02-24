#!/usr/bin/env python3
"""
S1_G: T Cell Subtype Marker Dotplot
Shows canonical markers for CD4+ T, CD8+ T, and NK cell subtypes.
Style matched to S1_F (flat gene list, clean layout).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import FULL_DATASET_H5AD

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib

# =============================================================================
# Configuration
# =============================================================================
SCALE = 4
DPI = 300
OUTPUT_DIR = Path(__file__).parent
OUTPUT_PNG = OUTPUT_DIR / "panel_S1_G.png"

# T/NK-specific markers (dict kept for reference, flattened for plotting)
markers = {
    'Pan-T': ['CD3D', 'CD3E'],
    'CD4': ['CD4', 'IL7R'],
    'CD8': ['CD8A', 'CD8B'],
    'NK': ['NKG7', 'GNLY', 'NCAM1'],
    'Treg': ['FOXP3', 'IL2RA'],
    'Cytotoxic': ['GZMB', 'PRF1'],
    'Exhaustion': ['PDCD1', 'HAVCR2', 'LAG3'],
}

cell_type_order = ['CD4+ T cells', 'CD8+ T cells', 'NK cells']

# =============================================================================
# Load data
# =============================================================================
print("Loading data...")
adata = sc.read_h5ad(FULL_DATASET_H5AD)
adata = adata[adata.obs['major_cell_type'].isin(cell_type_order)].copy()
print(f"Filtered to {adata.n_obs} T/NK cells")

# Verify markers
all_genes = [g for genes in markers.values() for g in genes]
available = [g for g in all_genes if g in adata.var_names]
missing = [g for g in all_genes if g not in adata.var_names]
if missing:
    print(f"Missing genes: {missing}")
    markers = {k: [g for g in v if g in adata.var_names] for k, v in markers.items()}
    markers = {k: v for k, v in markers.items() if v}
print(f"Using {len(available)}/{len(all_genes)} markers")

# =============================================================================
# Create dotplot — S1_F style (flat gene list, no return_fig)
# =============================================================================
plt.rcParams.update({
    'font.size': 4 * SCALE,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

fig_width = max(14, len(available) * 0.7)
fig_height = max(6, len(cell_type_order) * 1.2)
plt.figure(figsize=(fig_width, fig_height))

sc.pl.dotplot(
    adata,
    var_names=available,
    groupby='major_cell_type',
    categories_order=cell_type_order,
    use_raw=(adata.raw is not None),
    cmap='Reds',
    show=False,
    save=None,
    standard_scale='var',
)

plt.tight_layout()

print(f"Saving to {OUTPUT_PNG}")
plt.savefig(OUTPUT_PNG, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.savefig(OUTPUT_PNG.with_suffix('.svg'), dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close('all')
print("Done!")
