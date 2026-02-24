#!/usr/bin/env python3
"""
Panel F: MoMac Density Plots (2×2 grid)
Shows cell density for Pre R, Pre NR, Post R, Post NR conditions
4× scaling method - larger plots, smaller legend
"""

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.stats import gaussian_kde
from pathlib import Path
import sys

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *

# 4× scaling method
SCALE = 4
DPI = 300

# Target print size: 5.8 × 6.4 cm (fits Row 3 height)
FIGURE_WIDTH_CM = 5.8 * SCALE
FIGURE_HEIGHT_CM = 6.4 * SCALE

# Font sizes
TITLE_FONTSIZE = 5 * SCALE
LEGEND_FONTSIZE = 4 * SCALE  # Smaller legend
TICK_FONTSIZE = 4 * SCALE

# Point size (larger for visibility)
POINT_SIZE = 1.0 * SCALE

CMAP_DENSITY = 'magma'

# Condition definitions for 2x2 grid
CONDITIONS = [
    {'name': 'Pre R', 'phase': 'Pre', 'col': 'stomach_pre_grouping', 'val': 'Responsed'},
    {'name': 'Pre NR', 'phase': 'Pre', 'col': 'stomach_pre_grouping', 'val': 'No-response'},
    {'name': 'Post R', 'phase': 'Post', 'col': 'stomach_post_grouping', 'val': 'Responsed'},
    {'name': 'Post NR', 'phase': 'Post', 'col': 'stomach_post_grouping', 'val': 'No-response'},
]

def compute_density(x, y):
    """Compute KDE density for points."""
    xy = np.vstack([x, y])
    try:
        kde = gaussian_kde(xy)
        density = kde(xy)
    except:
        density = np.ones(len(x))
    return density

def main():
    print("=" * 60)
    print("Panel F: MoMac Density Plots (4× scaling)")
    print("=" * 60)

    script_dir = Path(__file__).parent
    data_path = MOMAC_H5AD

    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': TITLE_FONTSIZE,
        'axes.linewidth': 0.8 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    print("\nLoading data...")
    adata = sc.read_h5ad(data_path)
    print(f"  Total cells: {adata.n_obs:,}")

    # Get global UMAP limits
    umap = adata.obsm['X_umap']
    x_min, x_max = umap[:, 0].min(), umap[:, 0].max()
    y_min, y_max = umap[:, 1].min(), umap[:, 1].max()
    x_pad = (x_max - x_min) * 0.05
    y_pad = (y_max - y_min) * 0.05

    # Create figure
    fig_w = FIGURE_WIDTH_CM / 2.54
    fig_h = FIGURE_HEIGHT_CM / 2.54
    fig = plt.figure(figsize=(fig_w, fig_h))

    # GridSpec: 2×2 grid with small colorbar row
    gs = gridspec.GridSpec(3, 2, figure=fig,
                           height_ratios=[1, 1, 0.08],
                           hspace=0.12 * SCALE, wspace=0.08 * SCALE)

    print("\nCreating density plots...")
    scatter_handles = []

    for idx, cond in enumerate(CONDITIONS):
        row, col = idx // 2, idx % 2
        ax = fig.add_subplot(gs[row, col])

        # Filter cells
        mask = (
            (adata.obs['Treatment phase'] == cond['phase']) &
            (adata.obs[cond['col']] == cond['val'])
        )
        adata_sub = adata[mask]
        print(f"  {cond['name']}: {adata_sub.shape[0]:,} cells")

        x = adata_sub.obsm['X_umap'][:, 0]
        y = adata_sub.obsm['X_umap'][:, 1]

        # Compute density
        density = compute_density(x, y)
        order = np.argsort(density)
        x, y, density = x[order], y[order], density[order]

        # Plot density
        sc_plot = ax.scatter(x, y, c=density, s=POINT_SIZE,
                             cmap=CMAP_DENSITY, rasterized=False, linewidths=0)
        scatter_handles.append(sc_plot)

        ax.set_xlim(x_min - x_pad, x_max + x_pad)
        ax.set_ylim(y_min - y_pad, y_max + y_pad)
        ax.set_title(cond['name'], fontsize=TITLE_FONTSIZE, pad=2 * SCALE)
        ax.set_aspect('equal')
        ax.axis('off')

    # Colorbar (smaller, horizontal at bottom)
    ax_cbar = fig.add_subplot(gs[2, :])
    ax_cbar.axis('off')
    cbar = fig.colorbar(scatter_handles[0], ax=ax_cbar, orientation='horizontal',
                        fraction=0.6, pad=0.05, aspect=25)
    cbar.set_label('Cell Density', fontsize=LEGEND_FONTSIZE)
    cbar.ax.tick_params(labelsize=TICK_FONTSIZE * 0.8)

    plt.tight_layout()

    # Save
    output_path = script_dir / 'momac_density_panel.png'
    fig.savefig(output_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    fig.savefig(output_path.with_suffix('.svg'), bbox_inches='tight', facecolor='white')
    print(f"\nSaved: {output_path}")

    plt.close()

if __name__ == '__main__':
    main()
