#!/usr/bin/env python3
"""
Panel A: Large square UMAP of epithelial cells colored by minor cell states.
4× scaling method. Uses sc.pl.umap with right-margin legend (Figure 1 style).
"""

import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
from adjustText import adjust_text
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *

DPI = 300
SCALE = 4
CM_TO_INCH = 1 / 2.54

PANEL_SIZE_CM = 7.0 * SCALE

BASE_DIR = Path(__file__).parent

# Set2 palette (ColorBrewer) — ordered by category
COLORS = {
    'C0_Epi_PTMA': '#66C2A5',
    'C1_Epi_KRT19': '#FC8D62',
    'C2_Epi_CEACAM6': '#8DA0CB',
    'C3_Epi_Chief_Like_PGC': '#E78AC3',
    'C4_Epi_MUC5AC': '#A6D854',
    'C5_Epi_Stem_Like_TPX2': '#FFD92F',
    'C6_Epi_Stem_Like_SPINK4': '#E5C494',
    'C7_Epi_MT1E': '#B3B3B3',
    'C8_Epi_CD74': '#8DD3C7',
}

SHORT_LABELS = {
    'C0_Epi_PTMA': 'PTMA',
    'C1_Epi_KRT19': 'KRT19',
    'C2_Epi_CEACAM6': 'CEACAM5/6',
    'C3_Epi_Chief_Like_PGC': 'Chief_Like',
    'C4_Epi_MUC5AC': 'MUC5AC',
    'C5_Epi_Stem_Like_TPX2': 'Stem_TPX2',
    'C6_Epi_Stem_Like_SPINK4': 'Stem_SPINK4',
    'C7_Epi_MT1E': 'MT1E',
    'C8_Epi_CD74': 'CD74',
}


def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 5 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    print("Loading data...")
    adata = sc.read_h5ad(EPITHELIAL_H5AD)
    print(f"Loaded {adata.n_obs} cells")

    # Rename categories to short labels for legend
    orig_cats = list(adata.obs['minor_cell_state'].cat.categories)
    color_list = [COLORS[c] for c in orig_cats]
    adata.obs['minor_cell_state'] = adata.obs['minor_cell_state'].cat.rename_categories(SHORT_LABELS)
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

    # Add on-plot centroid labels with adjustText to prevent overlap
    coords = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata.obs_names)
    coords['cluster'] = adata.obs['minor_cell_state'].values

    texts = []
    for cluster_name in adata.obs['minor_cell_state'].cat.categories:
        mask = coords['cluster'] == cluster_name
        cx = coords.loc[mask, 'UMAP1'].median()
        cy = coords.loc[mask, 'UMAP2'].median()
        t = ax.text(cx, cy, cluster_name, fontsize=2.5 * SCALE, fontweight='bold',
                    ha='center', va='center',
                    bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                              alpha=0.8, edgecolor='none', linewidth=0))
        texts.append(t)

    # Auto-repel overlapping labels; draw thin leader lines from label to centroid
    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle='-', color='0.4', lw=0.4 * SCALE),
                expand=(1.5, 1.5),
                force_text=(0.8, 0.8),
                force_static=(0.5, 0.5),
                only_move={'text': 'xy'},
                ensure_inside_axes=True)

    # Rasterize scatter dots (keeps axes/legend as vectors, dots as embedded raster)
    for coll in ax.collections:
        coll.set_rasterized(True)

    output = BASE_DIR / 'epithelial_umap_minor_states.png'
    plt.savefig(output, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.savefig(output.with_suffix('.svg'), format='svg', dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {output}")


if __name__ == '__main__':
    main()
