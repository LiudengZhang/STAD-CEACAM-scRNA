#!/usr/bin/env python3
"""
Panel (right of Rows 1-2): UMAP of CD8+ T cells colored by minor cell states.
4× scaling method. Uses sc.pl.umap with right-margin legend (Figure 1 style).
"""

import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import TCD8_H5AD

SCALE = 4
DPI = 300
CM_TO_INCH = 1 / 2.54

PANEL_SIZE_CM = 7.0 * SCALE

BASE_DIR = Path(__file__).parent

# Set2 palette (ColorBrewer) — 8 CD8 states
COLORS = {
    'C0_CD8_Cytotoxic_CCL':   '#66C2A5',
    'C1_CD8_Cytotoxic_DUSP1': '#FC8D62',
    'C2_CD8_MAIT_KLRB1':      '#8DA0CB',
    'C3_CD8_Tcm_CCR7':        '#E78AC3',
    'C4_CD8_Temra_KLRG1':     '#A6D854',
    'C5_CD8_Prolif_MKI67':    '#FFD92F',
    'C6_CD8_Tex_PDCD1':       '#E5C494',
    'C7_CD8_ISG_ISG15':       '#B3B3B3',
}

SHORT_NAMES = {
    'C0_CD8_Cytotoxic_CCL':   'Cytotoxic_CCL',
    'C1_CD8_Cytotoxic_DUSP1': 'Cytotoxic_DUSP1',
    'C2_CD8_MAIT_KLRB1':      'MAIT',
    'C3_CD8_Tcm_CCR7':        'Tcm',
    'C4_CD8_Temra_KLRG1':     'Temra',
    'C5_CD8_Prolif_MKI67':    'Prolif',
    'C6_CD8_Tex_PDCD1':       'Tex',
    'C7_CD8_ISG_ISG15':       'ISG',
}


def main():
    print("=" * 60)
    print("CD8+ T cell UMAP (4× scaling, right-margin legend)")
    print("=" * 60)

    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 5 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    print("\nLoading data...")
    adata = sc.read_h5ad(TCD8_H5AD)
    print(f"  Total cells: {adata.n_obs:,}")

    # Rename categories to short labels for legend
    orig_cats = list(adata.obs['minor_cell_state'].cat.categories)
    color_list = [COLORS[c] for c in orig_cats]
    adata.obs['minor_cell_state'] = adata.obs['minor_cell_state'].cat.rename_categories(SHORT_NAMES)
    adata.uns['minor_cell_state_colors'] = color_list

    print("=== Color mapping ===")
    for cat, color in zip(adata.obs['minor_cell_state'].cat.categories, color_list):
        print(f"  {cat}: {color}")

    fig_size = PANEL_SIZE_CM * CM_TO_INCH
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))

    sc.pl.umap(
        adata,
        color='minor_cell_state',
        ax=ax,
        show=False,
        legend_loc='none',
        title='',
        frameon=False,
        size=8,
        alpha=0.6,
    )

    # Remove legend if scanpy created one anyway
    if ax.get_legend() is not None:
        ax.get_legend().remove()

    # Add on-plot centroid labels
    coords = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata.obs_names)
    coords['cluster'] = adata.obs['minor_cell_state'].values
    for cluster_name in adata.obs['minor_cell_state'].cat.categories:
        mask = coords['cluster'] == cluster_name
        cx = coords.loc[mask, 'UMAP1'].median()
        cy = coords.loc[mask, 'UMAP2'].median()
        ax.text(cx, cy, cluster_name, fontsize=3.5 * SCALE, fontweight='bold',
                ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.15', facecolor='white', alpha=0.7, edgecolor='none'))

    # Rasterize scatter dots (keeps axes/legend as vectors, dots as embedded raster)
    for coll in ax.collections:
        coll.set_rasterized(True)

    output_path = BASE_DIR / 'cd8_umap_minor_states.png'
    plt.savefig(output_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.savefig(output_path.with_suffix('.svg'), format='svg', dpi=DPI, bbox_inches='tight', facecolor='white')
    print(f"\nSaved: {output_path}")

    plt.close()


if __name__ == '__main__':
    main()
