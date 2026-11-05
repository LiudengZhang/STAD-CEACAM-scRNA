#!/usr/bin/env python3
"""
Figure 2A: epithelial UMAP coloured by minor cell state, with the nine
cluster names on the plot.

Reads the epithelial h5ad (EPITHELIAL_H5AD in 00_Config/paths.py).
Writes epithelial_umap_minor_states.[svg|pdf|png] beside this script.

Run: python create_epithelial_umap.py

Label positions and the short leader lines are fixed in data coordinates
(PINNED_LABELS, PINNED_LEADERS) rather than recomputed with adjustText, whose
result depends on the canvas size. A check at draw time stops the script if a
label is no longer nearest its own cluster.
"""

import numpy as np
import scanpy as sc
import pandas as pd
from pathlib import Path
from matplotlib.patches import FancyArrowPatch
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import EPITHELIAL_H5AD                          # noqa: E402
from cnsfig import Panel, style  # noqa: E402
from cnsfig.layout import umap_corner_axis  # noqa: E402


SCALE = 4                           # sets the point size relative to the tick size
SMALL_PT = 2.5

PANEL = Panel(2, "A", __file__, small_pt=SMALL_PT)
MARGIN = dict(left=2.9, right=0.6, top=0.6, bottom=0.6)

# Set2 palette (ColorBrewer) - ordered by category
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
    'C0_Epi_PTMA': 'Epi_PTMA',
    'C1_Epi_KRT19': 'Epi_KRT19',
    'C2_Epi_CEACAM6': 'Epi_CEACAM5/6',
    'C3_Epi_Chief_Like_PGC': 'Epi_Chief_Like',
    'C4_Epi_MUC5AC': 'Epi_MUC5AC',
    'C5_Epi_Stem_Like_TPX2': 'Epi_TPX2',
    'C6_Epi_Stem_Like_SPINK4': 'Epi_SPINK4',
    'C7_Epi_MT1E': 'Epi_MT1E',
    'C8_Epi_CD74': 'Epi_CD74',
}

# label -> position in UMAP data coordinates.
PINNED_LABELS = {
    'Epi_PTMA':        (-1.3673487686183599, 7.648700458808822),
    'Epi_KRT19':       (-0.6530729121408712, 2.522255165382308),
    'Epi_CEACAM5/6':   (1.5870579527498947, 1.1080370972594444),
    'Epi_Chief_Like':  (-4.344364528786972, -0.3045438816109474),
    'Epi_MUC5AC':      (-0.6054549738708204, -3.412629859642106),
    'Epi_TPX2':        (-2.0646697236417086, 11.30056164960472),
    'Epi_SPINK4':      (5.021559543478653, 3.3731412718733953),
    'Epi_MT1E':        (0.4446967487840592, 4.668013317390365),
    'Epi_CD74':        (2.840957321673275, 8.887738926216045),
}

# label -> (start, end) of its leader line, in UMAP data coordinates.
PINNED_LEADERS = {
    'Epi_PTMA':    ((-0.8578974478527783, 7.901319242287869),
                    (-0.851193573202325, 7.904643455200397)),
    'Epi_KRT19':   ((-0.0970986603142503, 2.7789300599185704),
                    (-0.09501596318800054, 2.779891572189104)),
    'Epi_MUC5AC':  ((0.11218166749815417, -3.152955187494662),
                    (0.12065659524927241, -3.1498885602373825)),
    'Epi_SPINK4':  ((5.3880739796948145, 3.6500102191982258),
                    (5.406764984130859, 3.6564266681671125)),
    'Epi_MT1E':    ((0.9270839981257755, 4.92634788467403),
                    (0.9736745953559858, 4.951298713684082)),
    'Epi_CD74':    ((2.608670181094002, 9.14865136248212),
                    (2.648776054382324, 9.171024322509766)),
}

LEADER_TOL = 0.15                   # data units


def _guard_pins(centroids):
    """Raise unless the pins still describe the data that was just loaded."""
    names = list(centroids)
    for name, (lx, ly) in PINNED_LABELS.items():
        d = {n: float(np.hypot(lx - centroids[n][0], ly - centroids[n][1]))
             for n in names}
        nearest = min(d, key=d.get)
        if nearest != name:
            raise SystemExit(
                f"pinned label {name!r} at ({lx:.4f}, {ly:.4f}) is nearer the "
                f"{nearest!r} cluster ({d[nearest]:.4f}) than its own "
                f"({d[name]:.4f}). The pins no longer describe this data.")
    for name, (_, end) in PINNED_LEADERS.items():
        cx, cy = centroids[name]
        d = float(np.hypot(end[0] - cx, end[1] - cy))
        if d > LEADER_TOL:
            raise SystemExit(
                f"leader for {name!r} ends {d:.4f} data units from that "
                f"cluster's median, over the {LEADER_TOL} tolerance.")


def main():
    style.apply(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2

    adata = sc.read_h5ad(EPITHELIAL_H5AD)

    # Rename categories to short labels for legend
    orig_cats = list(adata.obs['minor_cell_state'].cat.categories)
    color_list = [COLORS[c] for c in orig_cats]
    adata.obs['minor_cell_state'] = adata.obs['minor_cell_state'].cat.rename_categories(SHORT_LABELS)
    adata.uns['minor_cell_state_colors'] = color_list

    fig, ax = style.subplots_mm(PANEL.w_mm, PANEL.h_mm)
    style.margins_mm(fig, **MARGIN)

    sc.pl.umap(
        adata,
        color='minor_cell_state',
        ax=ax,
        show=False,
        legend_loc='none',
        title='',
        frameon=False,
        # Small, translucent dots so the clusters do not overplot.
        size=1 * AREA,
        alpha=0.55,
    )

    # Remove legend if scanpy created one anyway
    if ax.get_legend() is not None:
        ax.get_legend().remove()


    coords = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata.obs_names)
    coords['cluster'] = adata.obs['minor_cell_state'].values

    centroids = {}
    for cluster_name in adata.obs['minor_cell_state'].cat.categories:
        mask = coords['cluster'] == cluster_name
        centroids[cluster_name] = (float(coords.loc[mask, 'UMAP1'].median()),
                                   float(coords.loc[mask, 'UMAP2'].median()))
    _guard_pins(centroids)

    for cluster_name in adata.obs['minor_cell_state'].cat.categories:
        cx, cy = PINNED_LABELS[cluster_name]
        ax.text(cx, cy, cluster_name, **style.umap_label_kw())
    # White under-layer so the names read on top of the cloud.
    style.umap_label_underlays(ax)

    # The leader lines, in the same order the labels are drawn in.
    for cluster_name in adata.obs['minor_cell_state'].cat.categories:
        if cluster_name not in PINNED_LEADERS:
            continue
        posA, posB = PINNED_LEADERS[cluster_name]
        ax.add_patch(FancyArrowPatch(
            posA=posA, posB=posB,
            arrowstyle='-', connectionstyle='arc3,rad=0.0',
            shrinkA=0, shrinkB=0,
            color='0.4', lw=style.RULE_PT,
            transform=ax.transData))

    # Dots as an embedded raster; the type and the frame stay vectors.
    for coll in ax.collections:
        coll.set_rasterized(True)

    # UMAP1/UMAP2 arrows at the lower left, drawn after the leaders.
    umap_corner_axis(ax)

    PANEL.finish(fig, 'epithelial_umap_minor_states')


if __name__ == '__main__':
    main()
