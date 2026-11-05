#!/usr/bin/env python3
"""Figure 4D: UMAP of the monocyte-macrophage compartment, coloured by minor
cell state, with each state named at its cluster median.

Input:  MOMAC_H5AD (00_Config/paths.py)
Output: momac_umap_minor_states.{svg,pdf,png} in this directory
Run:    python create_momac_umap.py
"""

import scanpy as sc
import matplotlib
import matplotlib.transforms as mtransforms
matplotlib.use('Agg')
import pandas as pd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *                       # noqa: E402,F401,F403
from cnsfig import Panel, style  # noqa: E402
from cnsfig.layout import umap_corner_axis  # noqa: E402

SCALE = 4                       # marker-size scaling constants
SMALL_PT = 5.0

PANEL = Panel(4, "D", __file__, small_pt=SMALL_PT)


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

#: (dx, dy) in points from each cluster's median for the on-plot names; only
#: names that would otherwise touch are moved.
LABEL_NUDGE_PT = {
    'MoMac_IL1B': (0, 0),
    'Mac_TREM2': (4, 0),        # clear of MoMac_IL1B above it
    'MoMac_Inter': (0, 0),
    'Mono_CD14': (0, -2),       # clear of MoMac_Inter above it
    'Mac_MT1G': (0, 0),
    'Mac_Prolif': (0, 0),
    'Mono_CD16': (0, 0),
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
    PANEL.apply()
    AREA = (style.tick_pt() / (SMALL_PT * SCALE)) ** 2

    adata = sc.read_h5ad(MOMAC_H5AD)

    orig_cats = list(adata.obs['minor_cell_state'].cat.categories)
    color_list = [COLORS[c] for c in orig_cats]
    adata.obs['minor_cell_state'] = adata.obs['minor_cell_state'].cat.rename_categories(SHORT_NAMES)
    adata.uns['minor_cell_state_colors'] = color_list

    fig, ax = style.subplots_mm(PANEL.w_mm, PANEL.h_mm)

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

    if ax.get_legend() is not None:
        ax.get_legend().remove()

    # UMAP1/UMAP2 arrows at the lower left, shared by all UMAP panels.
    umap_corner_axis(ax)

    coords = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata.obs_names)
    coords['cluster'] = adata.obs['minor_cell_state'].values
    for cluster_name in adata.obs['minor_cell_state'].cat.categories:
        mask = coords['cluster'] == cluster_name
        cx = coords.loc[mask, 'UMAP1'].median()
        cy = coords.loc[mask, 'UMAP2'].median()
        dx, dy = LABEL_NUDGE_PT.get(cluster_name, (0, 0))
        tr = mtransforms.offset_copy(ax.transData, fig=fig, x=dx, y=dy,
                                     units='points')
        ax.text(cx, cy, cluster_name, **style.umap_label_kw(transform=tr))
    # White under-layer behind the labels.
    style.umap_label_underlays(ax)

    # Rasterize the dots only; text and axes stay vector.
    for coll in ax.collections:
        coll.set_rasterized(True)

    ax.set_aspect('equal', adjustable='datalim')
    PANEL.fit(fig)
    PANEL.finish(fig, 'momac_umap_minor_states')


if __name__ == '__main__':
    main()
