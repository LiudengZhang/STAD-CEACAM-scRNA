#!/usr/bin/env python3
"""
Figure 3C: UMAP of the CD8+ T cells, coloured by minor cell state, with the
state names on the plot.

Reads the CD8 T-cell h5ad (TCD8_H5AD in 00_Config/paths.py).
Writes cd8_umap_minor_states.[svg|pdf|png] beside this script.

Run: python create_cd8_umap.py

Each name is set at its state's median UMAP position, shifted by
LABEL_NUDGE_MM where neighbouring labels would collide. Points are
rasterised; text stays vector.
"""

import scanpy as sc
import matplotlib.transforms as mtransforms
import pandas as pd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import TCD8_H5AD
from cnsfig import Panel, style  # noqa: E402
from cnsfig.layout import umap_corner_axis  # noqa: E402

PANEL = Panel(3, "C", __file__, small_pt=3.5)
MARK, AREA = PANEL.mark, PANEL.area

# Set2 palette (ColorBrewer) - 8 CD8 states
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

# Label offsets in mm, (dx, dy) with dy up the page, for labels whose
# median positions would collide. The medians themselves are unchanged.
LABEL_NUDGE_MM = {
    'TCD8_ISG':   (7.5, -1.0),
    'TCD8_Tex':   (4.5, -3.0),
    'TCD8_MAIT':  (-3.0, 2.3),   # clear of DUSP1
    'TCD8_DUSP1': (0.5,  1.5),
    'TCD8_Tcm':   (-4.2, -1.5),  # clear of CCL
    'TCD8_CCL':   (-1.0, -2.6),  # clear of Tcm
}

SHORT_NAMES = {
    'C0_CD8_Cytotoxic_CCL':   'TCD8_CCL',
    'C1_CD8_Cytotoxic_DUSP1': 'TCD8_DUSP1',
    'C2_CD8_MAIT_KLRB1':      'TCD8_MAIT',
    'C3_CD8_Tcm_CCR7':        'TCD8_Tcm',
    'C4_CD8_Temra_KLRG1':     'TCD8_Temra',
    'C5_CD8_Prolif_MKI67':    'TCD8_Prolif',
    'C6_CD8_Tex_PDCD1':       'TCD8_Tex',
    'C7_CD8_ISG_ISG15':       'TCD8_ISG',
}


def main():
    adata = sc.read_h5ad(TCD8_H5AD)

    # Rename categories to short labels for legend
    orig_cats = list(adata.obs['minor_cell_state'].cat.categories)
    color_list = [COLORS[c] for c in orig_cats]
    adata.obs['minor_cell_state'] = adata.obs['minor_cell_state'].cat.rename_categories(SHORT_NAMES)
    adata.uns['minor_cell_state_colors'] = color_list

    fig, ax = PANEL.subplots()

    sc.pl.umap(
        adata,
        color='minor_cell_state',
        ax=ax,
        show=False,
        legend_loc='none',
        title='',
        frameon=False,
        # Small, translucent dots so the clusters do not overplot.
        size=2 * AREA,
        alpha=0.55,
    )

    # Remove legend if scanpy created one anyway
    if ax.get_legend() is not None:
        ax.get_legend().remove()

    # UMAP1/UMAP2 arrows at the lower left.
    umap_corner_axis(ax)

    # On-plot centroid labels, at each state's median embedding position.
    coords = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata.obs_names)
    coords['cluster'] = adata.obs['minor_cell_state'].values
    for cluster_name in adata.obs['minor_cell_state'].cat.categories:
        mask = coords['cluster'] == cluster_name
        cx = coords.loc[mask, 'UMAP1'].median()
        cy = coords.loc[mask, 'UMAP2'].median()
        dx, dy = LABEL_NUDGE_MM.get(cluster_name, (0.0, 0.0))
        tr = mtransforms.offset_copy(ax.transData, fig=fig,
                                     x=dx / 25.4, y=dy / 25.4, units='inches')
        ax.text(cx, cy, cluster_name, **style.umap_label_kw(transform=tr))
    # White under-layer so the names read on top of the cloud.
    style.umap_label_underlays(ax)

    # Dots as an embedded raster; the type and the labels stay vectors.
    for coll in ax.collections:
        coll.set_rasterized(True)

    PANEL.fit(fig)
    PANEL.finish(fig, 'cd8_umap_minor_states')


if __name__ == '__main__':
    main()
