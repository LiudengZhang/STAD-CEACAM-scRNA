#!/usr/bin/env python3
"""Figure 4F: Gaussian KDE cell density over the MoMac UMAP in the four
treatment-and-response groups. Each group is coloured over its own density
range; the colour bar is taken from the first group (Pre-R).

Input:  MOMAC_H5AD (00_Config/paths.py)
Output: momac_density_panel.{svg,pdf,png} in this directory
Run:    python create_density_panel.py
"""

import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.stats import gaussian_kde
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *                       # noqa: E402,F401,F403
from cnsfig import Panel, style  # noqa: E402

SCALE = 4                       # marker-size scaling constants
SMALL_PT = 4.0

PANEL = Panel(4, "F", __file__, small_pt=SMALL_PT)

CMAP_DENSITY = 'magma'

CONDITIONS = [
    {'name': 'Pre-R', 'phase': 'Pre', 'col': 'stomach_pre_grouping', 'val': 'Responsed'},
    {'name': 'Pre-NR', 'phase': 'Pre', 'col': 'stomach_pre_grouping', 'val': 'No-response'},
    {'name': 'Post-R', 'phase': 'Post', 'col': 'stomach_post_grouping', 'val': 'Responsed'},
    {'name': 'Post-NR', 'phase': 'Post', 'col': 'stomach_post_grouping', 'val': 'No-response'},
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
    PANEL.apply()
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2

    adata = sc.read_h5ad(MOMAC_H5AD)

    umap = adata.obsm['X_umap']
    x_min, x_max = umap[:, 0].min(), umap[:, 0].max()
    y_min, y_max = umap[:, 1].min(), umap[:, 1].max()
    x_pad = (x_max - x_min) * 0.05
    y_pad = (y_max - y_min) * 0.05

    fig = style.figure_mm(PANEL.w_mm, PANEL.h_mm)

    # 2 x 2 maps above a thin colour-bar row
    gs = gridspec.GridSpec(3, 2, figure=fig,
                           height_ratios=[1, 1, 0.08],
                           hspace=0.18, wspace=0.04)

    scatter_handles = []

    for idx, cond in enumerate(CONDITIONS):
        row, col = idx // 2, idx % 2
        ax = fig.add_subplot(gs[row, col])

        mask = (
            (adata.obs['Treatment phase'] == cond['phase']) &
            (adata.obs[cond['col']] == cond['val'])
        )
        adata_sub = adata[mask]

        x = adata_sub.obsm['X_umap'][:, 0]
        y = adata_sub.obsm['X_umap'][:, 1]

        density = compute_density(x, y)
        order = np.argsort(density)
        x, y, density = x[order], y[order], density[order]

        sc_plot = ax.scatter(x, y, c=density, s=1.0 * SCALE * AREA,
                             cmap=CMAP_DENSITY, rasterized=False, linewidths=0)
        scatter_handles.append(sc_plot)

        ax.set_xlim(x_min - x_pad, x_max + x_pad)
        ax.set_ylim(y_min - y_pad, y_max + y_pad)
        ax.set_title(cond['name'], pad=2 * SCALE * MARK)
        ax.set_aspect('equal')
        ax.axis('off')

    ax_cbar = fig.add_subplot(gs[2, :])
    ax_cbar.axis('off')
    cbar = fig.colorbar(scatter_handles[0], ax=ax_cbar, orientation='horizontal',
                        fraction=0.6, pad=0.05, aspect=25)
    cbar.set_label('Cell Density')
    cbar.ax.tick_params(labelsize=style.tick_pt())

    PANEL.fit(fig)
    PANEL.finish(fig, 'momac_density_panel')


if __name__ == '__main__':
    main()
