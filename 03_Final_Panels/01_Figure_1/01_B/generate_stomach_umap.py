#!/usr/bin/env python3
"""
Generate UMAP colored by Major Cell Type - Stomach Samples Only (32 samples)
============================================================================
Creates a UMAP visualization of stomach samples only, colored by major cell type.
"""

import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import FULL_DATASET_H5AD

# Set2 + Set3 palette (12 major cell types, consistent across Fig 1B & 1C)
CELL_TYPE_COLORS = {
    'Epithelial cells': '#e5c494',
    'T/NK cells': '#66c2a5',
    'Monocytes/Macrophages': '#fc8d62',
    'Plasma cells': '#8da0cb',
    'B cells': '#ffd92f',
    'Endothelial cells': '#a6d854',
    'Fibroblasts': '#e78ac3',
    'Neutrophils': '#b3b3b3',
    'Mast cells': '#fb8072',
    'Pericytes': '#bebada',
    'Dendritic cells': '#80b1d3',
    'Hepatocytes': '#bc80bd',
}

sc.set_figure_params(dpi=150, frameon=False, figsize=(12, 10), facecolor='white')
plt.rcParams.update({'svg.fonttype': 'none', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans']})

# Paths
DATA_FILE = FULL_DATASET_H5AD
OUTPUT_DIR = Path(__file__).parent

print("=" * 80)
print("GENERATING UMAP - STOMACH SAMPLES ONLY (32 samples)")
print("=" * 80)

# Load data
print(f"\nLoading data from: {DATA_FILE}")
adata = sc.read_h5ad(DATA_FILE)
print(f"Loaded: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# Filter to stomach samples only
adata_stomach = adata[adata.obs['Sample site'] == 'Stomach'].copy()
print(f"Stomach samples: {adata_stomach.shape[0]:,} cells")
print(f"Number of samples: {adata_stomach.obs['sample'].nunique()}")

# Merge T/NK cells and standardize names
print("\nStandardizing cell type names...")
adata_stomach.obs['major_cell_type_merged'] = adata_stomach.obs['major_cell_type'].replace({
    'CD4+ T cells': 'T/NK cells',
    'CD8+ T cells': 'T/NK cells',
    'NK cells': 'T/NK cells',
    'DC cells': 'Dendritic cells',
    'Fibroblast': 'Fibroblasts',
    'Pericyte': 'Pericytes',
    'Hepatocyte': 'Hepatocytes'
})
print(f"Merged major cell types: {adata_stomach.obs['major_cell_type_merged'].nunique()}")

# Verify UMAP coordinates exist
if 'X_umap' not in adata_stomach.obsm:
    raise ValueError("UMAP coordinates not found in adata.obsm['X_umap']")

print(f"UMAP coordinates shape: {adata_stomach.obsm['X_umap'].shape}")

# Print cell type distribution
print("\nCell type distribution (stomach only):")
cell_counts = adata_stomach.obs['major_cell_type_merged'].value_counts().sort_values(ascending=False)
for cell_type, count in cell_counts.items():
    pct = 100 * count / adata_stomach.shape[0]
    print(f"  {cell_type:25s}: {count:7,d} cells ({pct:5.2f}%)")

# Set up colors
cell_type_order = list(cell_counts.index)
colors = [CELL_TYPE_COLORS.get(ct, '#999999') for ct in cell_type_order]
adata_stomach.uns['major_cell_type_merged_colors'] = colors

# Create UMAP plot
print("\nCreating UMAP visualization...")
fig, ax = plt.subplots(figsize=(12, 10))
sc.pl.umap(
    adata_stomach,
    color='major_cell_type_merged',
    ax=ax,
    show=False,
    legend_loc='none',
    title='',
    frameon=False,
    size=2,
    alpha=0.6
)

# Remove legend if scanpy created one anyway
if ax.get_legend() is not None:
    ax.get_legend().remove()

# Add on-plot centroid labels
DISPLAY_NAMES = {
    'Epithelial cells': 'Epithelial',
    'T/NK cells': 'T/NK',
    'Monocytes/Macrophages': 'Mono/Mac',
    'Plasma cells': 'Plasma',
    'B cells': 'B cells',
    'Endothelial cells': 'Endothelial',
    'Fibroblasts': 'Fibroblasts',
    'Neutrophils': 'Neutrophils',
    'Mast cells': 'Mast',
    'Pericytes': 'Pericytes',
    'Dendritic cells': 'DC',
    'Hepatocytes': 'Hepatocytes',
}
coords = pd.DataFrame(adata_stomach.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata_stomach.obs_names)
coords['cluster'] = adata_stomach.obs['major_cell_type_merged'].values
for cluster_name in coords['cluster'].unique():
    mask = coords['cluster'] == cluster_name
    cx = coords.loc[mask, 'UMAP1'].median()
    cy = coords.loc[mask, 'UMAP2'].median()
    label = DISPLAY_NAMES.get(cluster_name, cluster_name)
    ax.text(cx, cy, label, fontsize=8, fontweight='bold',
            ha='center', va='center',
            bbox=dict(boxstyle='round,pad=0.15', facecolor='white', alpha=0.7, edgecolor='none'))

plt.tight_layout()

# Save
output_file = OUTPUT_DIR / 'umap_major_cell_type_stomach.png'
print(f"\nSaving to: {output_file}")
plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(OUTPUT_DIR / 'umap_major_cell_type_stomach.svg', format='svg', bbox_inches='tight', facecolor='white')
plt.savefig(OUTPUT_DIR / 'umap_major_cell_type_stomach.pdf', format='pdf', bbox_inches='tight')
plt.close()

print(f"✓ Saved: {output_file}")
print("=" * 80)
