#!/usr/bin/env python3
"""
Figure 4 panel D, RESTYLED (Version B) - UMAP of the MoMac compartment coloured
by minor cell state, with the seven state labels drawn on the plot.

Note the letter. This script lives in the directory `04_E`, but PROVENANCE.csv
records it as printed panel **D**: the one directory 04_E holds two printed
panels, D (this UMAP) and E (`create_momac_marker_dotplot.py`), and there is no
directory 04_D. CLAUDE.md rule 2 - never infer a panel letter from a directory
name.

Version A is
`03_Final_Panels/04_Figure_4/04_E/create_momac_umap.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every cell, every embedding coordinate, every category, every colour and every
label string is Version A's. The drawing code is the same code.

  printed panel  Figure 4 D     (PROVENANCE.csv; NOT inferred from "04_E")
  printed rect   33.5 x 42.6 mm    (panel_rects.csv; the square drawing is
                                    letterboxed inside it and prints ~33.5 mm)
  Version B box  62.0 x 62.0 mm

MARK
    Version A drew a 7 x 7 cm panel at SCALE = 4 (28 x 28 cm of canvas) and set
    its on-plot state labels at `3.5 * SCALE`, its smallest body type. The
    assembler fitted the saved 629.5 x 625.5 pt SVG into 33.5 mm of width, a
    fit of 0.1509, so that 14 pt type printed at 2.11 pt. So SMALL_PT = 3.5,
    SCALE = 4 and

        MARK = tick_pt / (SMALL_PT * SCALE) = 7 / 14 = 0.5
        AREA = MARK ** 2 = 0.25

    The one length Version A set is the scanpy dot `size=4`, a point area, so
    it carries `* AREA`. The label's `boxstyle='round,pad=0.15'` is in em, not
    points, and carries over unchanged.

THE BOX
    The state labels are the constraint. 'MoMac_Inter' is about 15 mm wide at
    7 pt; holding the printed label-to-panel proportion would need a 110 mm
    UMAP, which the page cannot give. 62 mm is the smallest square in which the
    seven labels sit on their own centroids without colliding and nothing
    overflows the canvas - the size PANEL_SPEC anticipates for a UMAP with
    on-plot labels. The panel is square because Version A's drawing is square;
    the published rect is taller only because the assembler letterboxed it.

TYPE
    `fontweight='bold'` is dropped from the state labels: cnsplots bolds axis
    titles and panel letters and nothing else (PANEL_SPEC, stage 1 finding 2).

Input : MOMAC_H5AD (00_Config/paths.py)
Output: this directory / momac_umap_minor_states.{svg,pdf,png}
"""

import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import *                       # noqa: E402,F401,F403
import panel_style_cns as style           # noqa: E402

# Version A's canvas convention, kept only so MARK can be derived from it.
SCALE = 4
SMALL_PT = 3.5                            # Version A's on-plot label `3.5 * SCALE`

PRINTED_MM = (33.5, 42.6)                 # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 62.0, 62.0
# 1 mm all round: the UMAP carries no axis furniture (frameon=False), so the
# only ink is the scatter and the seven on-plot labels, and both are inside the
# data limits. See the note beside style.overflow_mm below.
MARGIN = dict(left=1.0, right=1.0, top=1.0, bottom=1.0)

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

SHORT_NAMES = {
    'C0_Mac_Classic_TREM2':           'Mac_TREM2',
    'C1_Mono_Classic_CD14':           'Mono_CD14',
    'C2_MoMac_Intermediate_HLA-DRA':  'MoMac_Inter',
    'C3_Mac_Inflam_IL1B':             'Mac_IL1B',
    'C4_Mono_Alternative_CD16':       'Mono_CD16',
    'C5_Mac_Prolif_MKI67':            'Mac_Prolif',
    'C6_Mac_Metallothionein_MT1G':    'Mac_MT1G',
}


def main():
    print("=" * 60)
    print("Panel D: MoMac UMAP (restyled, drawn 1:1)")
    print("=" * 60)

    family = style.apply()
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

    # scanpy writes 'UMAP1'/'UMAP2' axis labels and then turns the axis off, so
    # the two Text artists stay live but are never drawn - `ax.axison` is False,
    # and neither string appears in any output file, in Version A or here. They
    # still count as ink to `style.overflow_mm`, and because the axis is never
    # drawn their positions are never updated either: the check reports a flat
    # 2.54 mm on the left and bottom whatever the margins are (measured at
    # 1, 3.6, 8 and 14 mm - the number does not move). That floor would mask a
    # real clip of up to 2.5 mm, which is exactly what the check exists to
    # catch, so the never-drawn artists are marked invisible. Nothing drawn
    # changes: `ax.get_xlabel()` still returns 'UMAP1' and the content gate
    # passes unchanged.
    ax.xaxis.label.set_visible(False)
    ax.yaxis.label.set_visible(False)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, BASE_DIR / 'momac_umap_minor_states')
    print(f"\nSaved: {BASE_DIR / 'momac_umap_minor_states'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
