#!/usr/bin/env python3
"""
Figure 1 panel B, RESTYLED (Version B) - UMAP of the stomach samples coloured
by major cell type.

Version A is
`03_Revised_Panels/Main_Figures/01_Figure_1/01_B/generate_stomach_umap.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows.
Every value read, every filter, every merge of cell-type names, every colour
and every label string is Version A's. The drawing code is the same code.

  printed panel  Figure 1 B     (PROVENANCE.csv; NOT inferred from "01_B")
  printed rect   69.4 x 56.8 mm   (panel_rects.csv)
  Version B box  66.0 x 55.0 mm

MARK
    Version A drew a 12 x 10 inch canvas - 304.8 x 254.0 mm - and set its only
    type, the on-plot centroid labels, at 8 pt. The assembler then fitted the
    saved SVG (299.2 mm wide) into 69.4 mm, a fit of 0.232, so those labels
    printed at 1.82 pt: the smallest type in Figure 1 and among the smallest in
    the paper. Version B draws 1:1 and the labels are set by cnsplots at 7 pt.

        MARK = tick_pt / SMALL_PT = 7 / 8 = 0.875
        AREA = MARK ** 2          = 0.766

    The one non-type size is scanpy's `size=2`, a marker area in pt^2. It
    printed at 0.32 pt across and now prints at 1.24 pt across - the same
    3.84x the type grew by, so the panel keeps its proportions.

    `bbox` padding on the label boxes is in units of the font size, so it is
    type and is left alone.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import FULL_DATASET_H5AD                      # noqa: E402
import panel_style_cns as style                          # noqa: E402

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

PRINTED_MM = (69.4, 56.8)           # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 66.0, 55.0
# scanpy draws its UMAP1/UMAP2 axis arrows and their labels outside the axes
# box when frameon=False; at 0.5 mm they ran 2.71 mm off the left and bottom
# edges. Found by style.overflow_mm, not by looking at it.
MARGIN = dict(left=3.4, right=0.5, top=0.5, bottom=3.4)

SMALL_PT = 8.0                      # Version A's only body type

sc.set_figure_params(dpi=150, frameon=False, figsize=(12, 10), facecolor='white')
# After sc.set_figure_params, so cnsplots' type system wins. This replaces
# Version A's local rcParams block, which asked for Arial (not installed here)
# and fell through to DejaVu Sans.
family = style.apply()
MARK = style.tick_pt() / SMALL_PT
AREA = MARK ** 2
print(f"  type set in {family}; body {style.body_pt():g} pt, "
      f"ticks/labels {style.tick_pt():g} pt; MARK {MARK:.3f}")

DATA_FILE = FULL_DATASET_H5AD
OUTPUT_DIR = Path(__file__).parent

print("=" * 80)
print("GENERATING UMAP - STOMACH SAMPLES ONLY (32 samples)")
print("=" * 80)

print(f"\nLoading data from: {DATA_FILE}")
adata = sc.read_h5ad(DATA_FILE)
print(f"Loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

adata_stomach = adata[adata.obs['Sample site'] == 'Stomach'].copy()
print(f"Stomach samples: {adata_stomach.shape[0]:,} cells")
print(f"Number of samples: {adata_stomach.obs['sample'].nunique()}")

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

if 'X_umap' not in adata_stomach.obsm:
    raise ValueError("UMAP coordinates not found in adata.obsm['X_umap']")

print(f"UMAP coordinates shape: {adata_stomach.obsm['X_umap'].shape}")

print("\nCell type distribution (stomach only):")
cell_counts = adata_stomach.obs['major_cell_type_merged'].value_counts().sort_values(ascending=False)
for cell_type, count in cell_counts.items():
    pct = 100 * count / adata_stomach.shape[0]
    print(f"  {cell_type:25s}: {count:7,d} cells ({pct:5.2f}%)")

cell_type_order = list(cell_counts.index)
colors = [CELL_TYPE_COLORS.get(ct, '#999999') for ct in cell_type_order]
adata_stomach.uns['major_cell_type_merged_colors'] = colors

print("\nCreating UMAP visualization...")
fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)
sc.pl.umap(
    adata_stomach,
    color='major_cell_type_merged',
    ax=ax,
    show=False,
    legend_loc='none',
    title='',
    frameon=False,
    size=2 * AREA,
    alpha=0.6
)

if ax.get_legend() is not None:
    ax.get_legend().remove()

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
    ax.text(cx, cy, label, fontsize=style.tick_pt(), fontweight='bold',
            ha='center', va='center',
            bbox=dict(boxstyle='round,pad=0.15', facecolor='white', alpha=0.7, edgecolor='none'))

style.margins_mm(fig, **MARGIN)
over = style.overflow_mm(fig)
if max(over) > 0.05:
    print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
          f"{tuple(round(v, 2) for v in over)}")
style.save_panel(fig, OUTPUT_DIR / 'umap_major_cell_type_stomach')
print(f"\nSaved: {OUTPUT_DIR / 'umap_major_cell_type_stomach'}.[svg|pdf|png] "
      f"at {PANEL_W_MM} x {PANEL_H_MM} mm")
print("=" * 80)
