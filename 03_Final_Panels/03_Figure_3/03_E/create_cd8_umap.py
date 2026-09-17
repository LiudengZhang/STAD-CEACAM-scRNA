#!/usr/bin/env python3
"""
Figure 3 panel C - UMAP of the CD8+ T cells, coloured by minor cell state.

The panel is drawn at the millimetre rectangle it prints in, read from
03_Final_Panels/panel_rects.csv through 00_Config/slots.py, so the type size
set here is the type size printed. Margins are measured from the rendered ink
rather than typed, and the panel-letter corner is left clear for the assembler.

  printed panel  Figure 3 C       (PROVENANCE.csv; the directory is 03_E, and
                                   the letter was looked up, not inferred)

Every cell, every embedding coordinate, every category, every colour and every
label string is the earlier drawing's: the same sc.pl.umap call with the same
alpha, the same category renaming, the same median-centroid label placement and
the same rasterisation of the dots.

MARK
    The earlier drawing used a canvas four times the printed size. The panel
    has no axes, no title and no legend, so its only body type is the on-plot
    centroid labels, set at 3.5 * SCALE. MARK carries the non-type point sizes
    across to the 1:1 canvas:

        MARK = style.tick_pt() / (SMALL_PT * SCALE)

    The one non-type size is scanpy's `size`, which is the matplotlib scatter
    area in points squared, so it is multiplied by AREA.

THE STATES CARRY THE PUBLISHED TCD8_ PREFIX AGAIN
    The centroid labels are placed at each state's median embedding position,
    so two neighbouring states are labelled a few millimetres apart whatever
    the panel is. At 44.0 mm of panel height "Cytotoxic_CCL" set 14.5 mm and
    printed over the label of the state beside it, and the labels were cut to
    the marker gene alone - CCL, DUSP1, MAIT.

    That was a loss the page did not have to take. The published panel prints
    TCD8_MAIT, TCD8_DUSP1 and six more, and the prefix is what tells the reader
    these are CD8 T-cell states at all. Since 2026-09-11 the panel is 54.0 mm
    tall (panel_rects_v2.csv), the embedding is that much further apart, and
    the widest of the eight - TCD8_DUSP1 at 13.7 mm - is narrower than the
    "Cytotoxic_CCL" that would not fit.

TYPE
    The centroid labels are set at the tick and legend size: they are the
    panel's smallest body type, the number MARK is derived from, and they do a
    legend's job on the plot. Their bold weight is dropped, because cnsplots
    bolds axis titles and panel letters and nothing else.
"""

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import pandas as pd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import TCD8_H5AD
import panel_style_cns as style  # noqa: E402
import slots  # noqa: E402

SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 3.5                  # the earlier smallest body type, before * SCALE

MARK = style.tick_pt() / (SMALL_PT * SCALE)   # length multiplier
AREA = MARK ** 2                              # area multiplier

PANEL_LETTER = "C"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(3, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(3, PANEL_LETTER)

BASE_DIR = Path(__file__).parent

# Set2 palette (ColorBrewer) — 8 CD8 states
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

#: Millimetres to move a centroid label by, (dx, dy) with dy up the page.
#
#  Two of the eight states are small and wedged between larger ones, so their
#  median embedding position lands within a label's width of a neighbour's:
#  the panel gate convicted TCD8_CCL against TCD8_ISG at 2.26 mm2. The
#  published page moves the same two labels clear of their centroids and leans
#  them on their own cluster, which is what these offsets reproduce. The
#  centroid itself is unchanged, and so is every cell that went into it.
LABEL_NUDGE_MM = {
    # At body size (2026-09-14) six labels are nudged off their medians so
    # that no two touch; each stays on its own cluster.
    'TCD8_ISG':   (7.5, -1.0),
    'TCD8_Tex':   (4.5, -3.0),
    'TCD8_MAIT':  (-2.5, 1.5),
    'TCD8_DUSP1': (0.5,  1.5),
    'TCD8_Tcm':   (-3.5, -1.5),
    'TCD8_CCL':   (-1.0, -2.0),
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
    print("=" * 60)
    print("Figure 3 panel C: CD8+ T cell UMAP")
    print("=" * 60)

    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

    print("\nLoading data...")
    adata = sc.read_h5ad(TCD8_H5AD)
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
        # AREA is the 4x-canvas area multiplier. 8 * AREA gives a 1.2 pt dot,
        # which at this cell count overplots to a flat block of Set2 at full
        # strength - the published panel is pastel because its dots are small
        # enough to leave the paper showing between them. Measured against the
        # published crop rather than chosen: a quarter of the area is half the
        # diameter.
        size=2 * AREA,
        alpha=0.55,
    )

    # Remove legend if scanpy created one anyway
    if ax.get_legend() is not None:
        ax.get_legend().remove()

    # scanpy sets the axis labels to 'UMAP1' / 'UMAP2' and then turns the whole
    # axis off (frameon=False), so neither is ever drawn - but both keep
    # visible=True and a position outside the axes, and the overflow check duly
    # reports millimetres of ink off the canvas that no reader will ever see.
    # Hiding the two artists makes the check see what is actually drawn. The
    # label text is untouched, and so is the image.
    ax.xaxis.label.set_visible(False)
    ax.yaxis.label.set_visible(False)

    # Add on-plot centroid labels
    coords = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata.obs_names)
    coords['cluster'] = adata.obs['minor_cell_state'].values
    for cluster_name in adata.obs['minor_cell_state'].cat.categories:
        mask = coords['cluster'] == cluster_name
        cx = coords.loc[mask, 'UMAP1'].median()
        cy = coords.loc[mask, 'UMAP2'].median()
        # Bold, black, and sitting straight on the embedding, which is how
        # the published panel sets them. The rounded white plate that was here
        # read as a row of buttons rather than as cluster labels, and it is not
        # on the page.
        dx, dy = LABEL_NUDGE_MM.get(cluster_name, (0.0, 0.0))
        tr = mtransforms.offset_copy(ax.transData, fig=fig,
                                     x=dx / 25.4, y=dy / 25.4, units='inches')
        # Body size, as the page's cluster labels are set (2026-09-14).
        ax.text(cx, cy, cluster_name, fontsize=style.body_pt(),
                ha='center', va='center', fontweight='bold', color='black',
                transform=tr)

    # Rasterize scatter dots (keeps axes/legend as vectors, dots as embedded raster)
    for coll in ax.collections:
        coll.set_rasterized(True)

    style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_CELL)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    style.save_panel(fig, BASE_DIR / 'cd8_umap_minor_states')
    print(f"\nSaved: {BASE_DIR / 'cd8_umap_minor_states'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
