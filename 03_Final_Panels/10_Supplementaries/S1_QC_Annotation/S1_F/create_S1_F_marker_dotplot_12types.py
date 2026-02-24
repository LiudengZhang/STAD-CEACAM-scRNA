#!/usr/bin/env python3
"""
S1_F: Cell Type Marker Dotplot (12 types — T/NK merged)
Canonical marker gene expression across 12 major cell types.
CD4+ T, CD8+ T, and NK cells are combined into "T/NK cells".
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
PANEL_WIDTH_CM = 22 * SCALE
PANEL_HEIGHT_CM = 10 * SCALE
OUTPUT_DIR = Path(__file__).parent
OUTPUT_PNG = OUTPUT_DIR / "panel_S1_F.png"

# Marker genes grouped by category
markers = {
    'T/NK': ['CD3D', 'CD3E', 'NKG7', 'GNLY'],
    'B cells': ['CD79A', 'MS4A1'],
    'Plasma': ['JCHAIN', 'MZB1'],
    'DC': ['FCER1A', 'CLEC10A'],
    'Mono/Mac': ['CD14', 'LYZ'],
    'Neutrophils': ['FCGR3B', 'CSF3R'],
    'Mast cells': ['TPSAB1', 'KIT'],
    'Epithelial': ['EPCAM', 'KRT18'],
    'Hepatocyte': ['ALB', 'APOA1'],
    'Endothelial': ['PECAM1', 'VWF'],
    'Fibroblast': ['COL1A1', 'DCN'],
    'Pericyte': ['RGS5', 'ACTA2'],
}

# 12 cell types (T/NK merged)
cell_type_map = {
    'CD4+ T cells': 'T/NK cells',
    'CD8+ T cells': 'T/NK cells',
    'NK cells': 'T/NK cells',
    'B cells': 'B cells',
    'Plasma cells': 'Plasma cells',
    'DC cells': 'DC cells',
    'Monocytes/Macrophages': 'Monocytes/Macrophages',
    'Neutrophils': 'Neutrophils',
    'Mast cells': 'Mast cells',
    'Epithelial cells': 'Epithelial cells',
    'Hepatocyte': 'Hepatocyte',
    'Endothelial cells': 'Endothelial cells',
    'Fibroblast': 'Fibroblast',
    'Pericyte': 'Pericyte',
}

cell_type_order = [
    'T/NK cells',
    'B cells',
    'Plasma cells',
    'DC cells',
    'Monocytes/Macrophages',
    'Neutrophils',
    'Mast cells',
    'Epithelial cells',
    'Hepatocyte',
    'Endothelial cells',
    'Fibroblast',
    'Pericyte',
]

# =============================================================================
# Load and prepare data
# =============================================================================
print("Loading data...")
adata = sc.read_h5ad(FULL_DATASET_H5AD)
print(f"Loaded {adata.n_obs} cells")

# Map to 12 types
adata.obs['cell_type_12'] = adata.obs['major_cell_type'].map(cell_type_map)
adata = adata[adata.obs['cell_type_12'].notna()].copy()
adata.obs['cell_type_12'] = adata.obs['cell_type_12'].astype('category')

# Verify markers exist
all_genes = [g for genes in markers.values() for g in genes]
available = [g for g in all_genes if g in adata.var_names]
missing = [g for g in all_genes if g not in adata.var_names]
if missing:
    print(f"Missing genes: {missing}")
    markers = {k: [g for g in v if g in adata.var_names] for k, v in markers.items()}
    markers = {k: v for k, v in markers.items() if v}
print(f"Using {len(available)}/{len(all_genes)} markers")

# =============================================================================
# Create dotplot
# =============================================================================
plt.rcParams.update({
    'font.size': 6 * SCALE,
    'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

fig_width = max(14, len(available) * 0.35)
fig_height = max(8, len(cell_type_order) * 0.6)
plt.figure(figsize=(fig_width, fig_height))

sc.pl.dotplot(
    adata,
    var_names=available,
    groupby='cell_type_12',
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
