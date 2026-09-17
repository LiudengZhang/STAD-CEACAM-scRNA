#!/usr/bin/env python3
"""
Figure 4 panel F - cell density over the MoMac embedding in the four
treatment-and-response groups, on one shared density scale.

  printed panel  Figure 4 F       (PROVENANCE.csv - the directory is "04_F";
                                   do NOT read the directory as the letter)

The embedding, the four group definitions, the kernel density estimate, the
draw order, the shared axis limits, the colour map and the shared colour bar
are unchanged. Only the canvas and the type change: the panel is drawn at the
millimetre rectangle it prints in and set in the figure's one type system.

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the colour-bar label and its tick labels - at
    `4 * SCALE`, so SMALL_PT = 4 and

        MARK = style.tick_pt() / (SMALL_PT * SCALE)
        AREA = MARK ** 2

    The point size is an area and is scaled by AREA; the title pad is a length
    in points and is scaled by MARK. The grid's `hspace` and `wspace` are
    fractions of an axes, not lengths, so they are carried over unchanged and
    the four boxes sit at the same relative spacing as before.

Judgement call, stated plainly:
  - the colour-bar tick labels were set at 0.8 of the tick size. At the type
    spec that is 4.8 pt, below the floor this figure is set to, so the factor
    is dropped and they take the tick size from the system.

Input : MOMAC_H5AD (00_Config/paths.py)
Output: this directory / momac_density_panel.{svg,pdf,png}
"""

import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.stats import gaussian_kde
from pathlib import Path
import sys

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *                       # noqa: E402,F401,F403
import panel_style_cns as style           # noqa: E402
import slots                              # noqa: E402

SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 4.0                  # the earlier smallest body type, before * SCALE

PANEL_LETTER = "F"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(4, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(4, PANEL_LETTER)

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
    print("Panel F: MoMac Density Plots")
    print("=" * 60)

    script_dir = Path(__file__).parent
    data_path = MOMAC_H5AD

    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.4f}")

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
    fig = style.figure_mm(PANEL_W_MM, PANEL_H_MM)

    # GridSpec: 2×2 grid with small colorbar row
    gs = gridspec.GridSpec(3, 2, figure=fig,
                           height_ratios=[1, 1, 0.08],
                           # Closed up 2026-09-14 (evening): the four maps
                           # were 0.48/0.32 apart, and the row gave the
                           # width saved to D and G (build_grid_v2.X_NEW).
                           hspace=0.18, wspace=0.04)

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
        sc_plot = ax.scatter(x, y, c=density, s=1.0 * SCALE * AREA,
                             cmap=CMAP_DENSITY, rasterized=False, linewidths=0)
        scatter_handles.append(sc_plot)

        ax.set_xlim(x_min - x_pad, x_max + x_pad)
        ax.set_ylim(y_min - y_pad, y_max + y_pad)
        ax.set_title(cond['name'], pad=2 * SCALE * MARK)
        ax.set_aspect('equal')
        ax.axis('off')

    # Colorbar (smaller, horizontal at bottom)
    ax_cbar = fig.add_subplot(gs[2, :])
    ax_cbar.axis('off')
    cbar = fig.colorbar(scatter_handles[0], ax=ax_cbar, orientation='horizontal',
                        fraction=0.6, pad=0.05, aspect=25)
    cbar.set_label('Cell Density')
    cbar.ax.tick_params(labelsize=style.tick_pt())

    style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_CELL)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    stem = 'momac_density_panel'
    style.save_panel(fig, script_dir / stem)
    print(f"\nSaved: {stem}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
