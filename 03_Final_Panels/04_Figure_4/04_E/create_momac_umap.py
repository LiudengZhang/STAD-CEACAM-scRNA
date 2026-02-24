#!/usr/bin/env python3
"""
Panel E: UMAP of MoMac cells colored by minor cell states.
4× scaling method. Uses sc.pl.umap with right-margin legend (Figure 1 style).
"""

import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *

SCALE = 4
DPI = 300
CM_TO_INCH = 1 / 2.54

PANEL_SIZE_CM = 7.0 * SCALE

BASE_DIR = Path(__file__).parent

# Set2 palette (ColorBrewer) — 7 MoMac states
COLORS = {
    'C0_Mac_Classic_TREM2': '#66C2A5',
    'C1_Mono_Classic_CD14': '#FC8D62',
    'C2_MoMac_Intermediate_HLA-DRA': '#8DA0CB',
    'C3_Mac_Inflam_IL1B': '#E78AC3',
    'C4_Mono_Alternative_CD16': '#A6D854',
    'C5_Mac_Prolif_MKI67': '#FFD92F',
    'C6_Mac_Metallothionein_MT1G': '#E5C494',
}

SHORT_NAMES = {
    'C0_Mac_Classic_TREM2':           'Mac_TREM2',
    'C1_Mono_Classic_CD14':           'Mono_CD14',
    'C2_MoMac_Intermediate_HLA-DRA':  'MoMac_Inter',
    'C3_Mac_Inflam_IL1B':             'Mac_IL1B',
    'C4_Mono_Alternative_CD16':       'Mono_CD16',
    'C5_Mac_Prolif_MKI67':            'Mac_Prolif',
    'C6_Mac_Metallothionein_MT1G':    'Mac_MT1G',
}


def main():
    print("=" * 60)
    print("Panel E: MoMac UMAP (4× scaling, right-margin legend)")
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
    adata = sc.read_h5ad(MOMAC_H5AD)
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
        size=4,
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

    output_path = BASE_DIR / 'momac_umap_minor_states.png'
    plt.savefig(output_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.savefig(output_path.with_suffix('.svg'), format='svg', dpi=DPI, bbox_inches='tight', facecolor='white')
    print(f"\nSaved: {output_path}")

    plt.close()


if __name__ == '__main__':
    main()
