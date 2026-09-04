#!/usr/bin/env python3
"""
Figure 3 panel C, RESTYLED (Version B) - UMAP of CD8+ T cells coloured by minor
cell state.

NOTE THE PANEL LETTER. This directory is `03_E` but PROVENANCE.csv says it
holds printed panel **C**. CLAUDE.md rule 2: never infer a panel letter from a
directory name.

Version A is
`03_Revised_Panels/Main_Figures/03_Figure_3/03_E/create_cd8_umap.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every cell, every embedding coordinate, every category, every colour and every
label string is Version A's. The drawing code is the same code: the same
`sc.pl.umap` call with the same `alpha`, the same category renaming, the same
median-centroid label placement and the same rasterisation of the dots.

  printed panel  Figure 3 C     (PROVENANCE.csv; NOT inferred from "03_E")
  printed rect   45.1 x 44.0 mm   (panel_rects.csv)
  Version B box  62.0 x 62.0 mm

MARK
    Version A drew a 28.0 x 28.0 cm canvas (SCALE = 4). The panel has no axes
    (`frameon=False`), no title and no legend, so its only body type is the
    on-plot centroid labels at `fontsize=3.5 * SCALE`. SMALL_PT = 3.5 and, by
    PANEL_SPEC.md,

        MARK = tick_pt / (SMALL_PT * SCALE) = 7 / 14 = 0.5
        AREA = MARK ** 2                    = 0.25

    The one non-type size is scanpy's `size=8`, which is the matplotlib scatter
    area in pt^2, so it is multiplied by AREA.

TYPE
    The centroid labels are set at `style.tick_pt()`. They are the panel's
    smallest body type - the number MARK is derived from - and they are doing a
    legend's job on the plot, which is cnsplots' 7 pt slot. `fontweight='bold'`
    is dropped: cnsplots bolds axis titles and panel letters and nothing else
    (PANEL_SPEC.md, Type).

    The box grew from 45.1 mm because the longest label, `Cytotoxic_DUSP1`,
    sets ~13 mm at 7 pt; at 45 mm it would span a third of the embedding.
"""

import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import TCD8_H5AD
import panel_style_cns as style  # noqa: E402

# Version A's canvas multiplier. Version B draws 1:1, so SCALE survives only to
# reproduce the exact numbers Version A set for its non-type lengths.
SCALE = 4
SMALL_PT = 3.5                       # Version A's smallest body type

# See MARK above. `style.tick_pt()` reads cnsplots' own setting, so the factor
# is derived, never a literal.
MARK = style.tick_pt() / (SMALL_PT * SCALE)   # 0.5, length multiplier
AREA = MARK ** 2                              # 0.25, area multiplier

PRINTED_MM = (45.1, 44.0)            # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 62.0, 62.0
# Millimetres of paper. There is no axis furniture at all; the margin is only
# the room the outermost centroid labels need to stay on the canvas.
MARGIN = dict(left=5.0, right=5.0, top=2.0, bottom=2.0)

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

SHORT_NAMES = {
    'C0_CD8_Cytotoxic_CCL':   'Cytotoxic_CCL',
    'C1_CD8_Cytotoxic_DUSP1': 'Cytotoxic_DUSP1',
    'C2_CD8_MAIT_KLRB1':      'MAIT',
    'C3_CD8_Tcm_CCR7':        'Tcm',
    'C4_CD8_Temra_KLRG1':     'Temra',
    'C5_CD8_Prolif_MKI67':    'Prolif',
    'C6_CD8_Tex_PDCD1':       'Tex',
    'C7_CD8_ISG_ISG15':       'ISG',
}


def main():
    print("=" * 60)
    print("CD8+ T cell UMAP (Version B, drawn 1:1 at print size)")
    print("=" * 60)

    family = style.apply()
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
        size=8 * AREA,
        alpha=0.6,
    )

    # Remove legend if scanpy created one anyway
    if ax.get_legend() is not None:
        ax.get_legend().remove()

    # scanpy sets the axis labels to 'UMAP1' / 'UMAP2' and then turns the whole
    # axis off (`frameon=False`), so neither is ever drawn - but both keep
    # `visible=True` and a position outside the axes, and `style.overflow_mm`
    # duly reports 2.5 mm of ink off the canvas that no reader will ever see.
    # Hiding the two artists makes the check see what is actually drawn. The
    # label *text* is untouched - the gate reads `ax.get_xlabel()`, which still
    # returns 'UMAP1' - and so is the image.
    ax.xaxis.label.set_visible(False)
    ax.yaxis.label.set_visible(False)

    # Add on-plot centroid labels
    coords = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata.obs_names)
    coords['cluster'] = adata.obs['minor_cell_state'].values
    for cluster_name in adata.obs['minor_cell_state'].cat.categories:
        mask = coords['cluster'] == cluster_name
        cx = coords.loc[mask, 'UMAP1'].median()
        cy = coords.loc[mask, 'UMAP2'].median()
        ax.text(cx, cy, cluster_name, fontsize=style.tick_pt(),
                ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.15', facecolor='white', alpha=0.7, edgecolor='none'))

    # Rasterize scatter dots (keeps axes/legend as vectors, dots as embedded raster)
    for coll in ax.collections:
        coll.set_rasterized(True)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, BASE_DIR / 'cd8_umap_minor_states')
    print(f"\nSaved: {BASE_DIR / 'cd8_umap_minor_states'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
