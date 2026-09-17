#!/usr/bin/env python3
"""
Figure 4 panel D - UMAP of the monocyte-macrophage compartment, coloured by
minor cell state, with each state named at its own centroid.

  printed panel  Figure 4 D       (PROVENANCE.csv - the directory is "04_E",
                                   which also holds the drawing for panel E;
                                   do NOT read the directory as the letter)

The embedding, the seven states, the colour assignment, the short state names,
the centroid positions and the point size expression are unchanged. Only the
canvas and the type change: the panel is drawn at the millimetre rectangle it
prints in and set in the figure's one type system.

ONE STATE NAME IS RESTORED
    The page names the C3 state MoMac_IL1B; the script still carried the
    earlier Mac_IL1B. The cells, the centroid and the colour do not move.

MARK
    The earlier drawing used a canvas four times the printed size and set its
    base type from `font.size = 5 * SCALE`, so SMALL_PT = 5 and

        MARK = style.tick_pt() / (SMALL_PT * SCALE)
        AREA = MARK ** 2

    The scatter's point size is an area and is scaled by AREA.

Judgement call, stated plainly:
  - the centroid labels were set at `3.5 * SCALE`, which under MARK is 4.2 pt,
    below the floor this figure is set to. They take their size from the system
    instead, at tick_pt. Their weight is kept: these strings are drawn on top
    of the point cloud, and the weight and the white background box behind them
    are together what separates them from it.

Input : MOMAC_H5AD (00_Config/paths.py)
Output: this directory / momac_umap_minor_states.{svg,pdf,png}
"""

import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *                       # noqa: E402,F401,F403
import panel_style_cns as style           # noqa: E402
import slots                              # noqa: E402

SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 5.0                  # the earlier base type, before * SCALE

PANEL_LETTER = "D"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(4, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(4, PANEL_LETTER)

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

#: (dx, dy) in points from each cluster's median, for the on-plot names.
LABEL_NUDGE_PT = {
    'MoMac_IL1B': (-8, 6),
    'Mac_TREM2': (10, 2),
    'MoMac_Inter': (0, -1),
    'Mono_CD14': (-6, -6),
    'Mac_MT1G': (8, -2),
    'Mac_Prolif': (14, -5),
    'Mono_CD16': (-6, -6),
}

SHORT_NAMES = {
    'C0_Mac_Classic_TREM2':           'Mac_TREM2',
    'C1_Mono_Classic_CD14':           'Mono_CD14',
    'C2_MoMac_Intermediate_HLA-DRA':  'MoMac_Inter',
    'C3_Mac_Inflam_IL1B':             'MoMac_IL1B',
    'C4_Mono_Alternative_CD16':       'Mono_CD16',
    'C5_Mac_Prolif_MKI67':            'Mac_Prolif',
    'C6_Mac_Metallothionein_MT1G':    'Mac_MT1G',
}


def main():
    print("=" * 60)
    print("Panel D: MoMac UMAP")
    print("=" * 60)

    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.4f}")

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

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

    sc.pl.umap(
        adata,
        color='minor_cell_state',
        ax=ax,
        show=False,
        legend_loc='none',
        title='',
        frameon=False,
        size=4 * AREA,
        alpha=0.6,
    )

    # Remove legend if scanpy created one anyway
    if ax.get_legend() is not None:
        ax.get_legend().remove()

    # scanpy names the two embedding axes and then switches the axis off, so
    # neither name has ever been drawn - the saved panel carries no such
    # string. The Text artists survive the switch, though, and report an
    # extent, which the 1:1 overflow check reads as ink 2.3 mm off two edges.
    # Marking them invisible tells the measurement what the page already
    # shows; the labels themselves are left in place, unchanged.
    ax.xaxis.label.set_visible(False)
    ax.yaxis.label.set_visible(False)

    # Add on-plot centroid labels
    coords = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata.obs_names)
    coords['cluster'] = adata.obs['minor_cell_state'].values
    for cluster_name in adata.obs['minor_cell_state'].cat.categories:
        mask = coords['cluster'] == cluster_name
        cx = coords.loc[mask, 'UMAP1'].median()
        cy = coords.loc[mask, 'UMAP2'].median()
        # Body size, bold, black, no plate: how the published page and the
        # other UMAPs of this paper set their cluster labels (2026-09-14).
        # Nudged off the medoid in points (2026-09-14, evening): at the
        # map's true aspect two pairs of names overprinted. The anchor is the
        # cluster's median; only the text moves.
        dx, dy = LABEL_NUDGE_PT.get(cluster_name, (0, 0))
        ax.annotate(cluster_name, (cx, cy), xytext=(dx, dy),
                    textcoords='offset points', fontsize=style.body_pt(),
                    fontweight='bold', ha='center', va='center', color='black')

    # Rasterize scatter dots (keeps axes/legend as vectors, dots as embedded raster)
    for coll in ax.collections:
        coll.set_rasterized(True)

    # Equal aspect (2026-09-14, evening): in a 33.5 x 47 mm slot the map was
    # stretched tall ("D is squashed" - the author). The row is re-split so D
    # is 38.6 mm wide, and the map keeps its own proportions inside the box.
    ax.set_aspect('equal', adjustable='datalim')
    style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_CELL)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    stem = 'momac_umap_minor_states'
    style.save_panel(fig, BASE_DIR / stem)
    print(f"\nSaved: {stem}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
