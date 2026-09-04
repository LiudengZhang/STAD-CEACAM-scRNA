#!/usr/bin/env python3
"""
Figure 4 panel F, RESTYLED (Version B) - MoMac UMAP cell density in the four
treatment/response groups (Pre R, Pre NR, Post R, Post NR).

Version A is
`03_Final_Panels/04_Figure_4/04_F/create_density_panel.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every cell selection, every KDE, every axis limit and every string is Version
A's. The drawing code is the same code.

  printed panel  Figure 4 F     (PROVENANCE.csv; NOT inferred from "04_F")
  printed rect   50.7 x 47.4 mm    (panel_rects.csv)
  Version B box  62.0 x 56.0 mm

MARK
    Version A drew 5.8 x 6.4 cm at SCALE = 4 and its smallest body type was
    `LEGEND_FONTSIZE = TICK_FONTSIZE = 4 * SCALE` (the panel titles were
    `5 * SCALE`). The assembler fitted the saved 524.1 x 614.2 pt SVG into
    50.7 x 47.4 mm, a fit of 0.2188, so that 16 pt type printed at 3.50 pt.
    So SMALL_PT = 4, SCALE = 4 and

        MARK = tick_pt / (SMALL_PT * SCALE) = 7 / 16 = 0.4375
        AREA = MARK ** 2 = 0.1914

    The two lengths Version A set in points are the density point size (an
    area, so AREA) and the title pad; both are written exactly as Version A
    wrote them with the factor appended. `hspace` and `wspace` are fractions of
    the axes size, not points, so they carry over unchanged.

TYPE
    Version A set the colorbar tick labels at `TICK_FONTSIZE * 0.8`, i.e. a
    nominal 3.2 pt that printed at 2.80 pt, and the colorbar label at
    `LEGEND_FONTSIZE`. Both explicit sizes are deleted: at 1:1 the panel has the
    room, and cnsplots' 7 pt tick / 8 pt label is the whole point of the
    restyle. Nothing else about the colorbar changes.

THE FRAME IS cnsplots'
    Version A's `axes.linewidth` is the axes frame, which cnsplots declares and
    `describe()` reports; it is dropped so the whole figure carries one frame
    weight. All four density axes are `axis('off')` in any case.

Input : MOMAC_H5AD (00_Config/paths.py)
Output: this directory / momac_density_panel.{svg,pdf,png}
"""

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.stats import gaussian_kde
from pathlib import Path
import sys

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import *                       # noqa: E402,F401,F403
import panel_style_cns as style           # noqa: E402

# Version A's canvas convention, kept only so MARK can be derived from it.
SCALE = 4
SMALL_PT = 4.0                            # Version A's `TICK_FONTSIZE = 4 * SCALE`

PRINTED_MM = (50.7, 47.4)                 # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 62.0, 56.0
MARGIN = dict(left=3.4, right=3.4, top=3.5, bottom=7.5)

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
    print("Panel F: MoMac Density Plots")
    print("=" * 60)

    script_dir = Path(__file__).parent
    data_path = MOMAC_H5AD

    family = style.apply()
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
        sc_plot = ax.scatter(x, y, c=density, s=POINT_SIZE * AREA,
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

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, script_dir / 'momac_density_panel')
    print(f"\nSaved: {script_dir / 'momac_density_panel'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")

if __name__ == '__main__':
    main()
