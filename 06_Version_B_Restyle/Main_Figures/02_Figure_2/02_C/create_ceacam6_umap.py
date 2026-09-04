#!/usr/bin/env python3
"""
Figure 2, printed panel C, RESTYLED (Version B) - CEACAM6 expression UMAP on
epithelial cells.

Version A is
`03_Final_Panels/02_Figure_2/02_C/create_ceacam6_umap.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md
allows: the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every cell, the gene, the colour map and both colour limits (vmin 2, vmax 6)
are Version A's. Expression is read from `.raw`, exactly as Version A reads it.
The drawing code is the same code.

  printed panel  Figure 2 C     (PROVENANCE.csv; the directory letter happens
                                 to agree here - it was still looked up. Note
                                 Version A's own docstring calls it "Panel 2D",
                                 which is the pre-submission lettering and is
                                 wrong; PROVENANCE.csv is the authority.)
  printed rect   24.7 x 23.1 mm    (panel_rects.csv)
  Version B box  62.0 x 52.0 mm

MARK
    Version A drew at SCALE = 4 (3.5 x 3.5 cm x 4 = 140 x 140 mm) and its
    smallest body type is the colorbar label and its tick labels, both at
    `5 * SCALE`. So

        SCALE = 4, SMALL_PT = 5
        MARK  = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 20 = 0.35
        AREA  = MARK ** 2 = 0.1225

    The one non-type size is the UMAP dot area, `size = 3 * SCALE`, which takes
    AREA.

    Version A's colorbar `tick_params(width=1.0 * SCALE, length=4 * SCALE)` and
    `outline.set_linewidth(1.0 * SCALE)` are NOT carried over: a colorbar's
    ticks and outline are axes furniture, cnsplots has its own settings for
    them (axes.linewidth 0.5, tick width 0.6, length 2), and following the
    library rather than rescaling the old numbers is the standard-methods rule.
    Nothing they control is a plotted value.

    The panel grew from 24.7 x 23.1 mm to 62 x 52 mm. The title and the
    colorbar legend at 7/8 pt need about 14 mm of the width before any cell is
    drawn, and 149,373 dots at their MARK-held size are less crowded in the
    larger square, not more.

LAYOUT ORDER
    `style.margins_mm` is called before `sc.pl.umap`, not after. scanpy builds
    the colorbar with `plt.colorbar(..., ax=ax)`, which steals its space out of
    the axes' *current* position; calling subplots_adjust afterwards would move
    the axes back over the colorbar. Version A got away with `tight_layout()`
    last because tight_layout knows about the colorbar's axes. Nothing drawn
    moves; only the order of two layout calls.
"""

import scanpy as sc
import numpy as np
import os

# Central config
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import *                                       # noqa: E402,F403
import panel_style_cns as style                           # noqa: E402

SCALE = 4                           # Version A's canvas multiplier
SMALL_PT = 5.0                      # Version A's smallest body type

PRINTED_MM = (24.7, 23.1)           # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 62.0, 52.0
MARGIN = dict(left=1.5, right=13.0, top=5.5, bottom=1.5)

# Paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = BASE_DIR


def main():
    family = style.apply()
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.5f}")

    print("Loading data...")
    adata = sc.read_h5ad(EPITHELIAL_H5AD)                 # noqa: F405
    print(f"Loaded {adata.n_obs} cells")

    gene = 'CEACAM6'
    # The working inputs keep the log1p matrix in .raw; the clean deposit
    # promotes it to .X and carries no .raw, so a file without .raw already
    # holds the same numbers in .X.
    src = adata.raw if adata.raw is not None else adata
    idx = src.var_names.get_loc(gene)
    expr = src.X[:, idx]
    if hasattr(expr, 'toarray'):
        expr = expr.toarray().flatten()
    else:
        expr = np.array(expr).flatten()
    adata.obs['CEACAM6_expr'] = expr

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)
    # Before the colorbar - see LAYOUT ORDER in the header.
    style.margins_mm(fig, **MARGIN)

    sc.pl.umap(
        adata,
        color='CEACAM6_expr',
        ax=ax,
        show=False,
        frameon=False,
        title='',
        size=3 * SCALE * AREA,
        cmap='Reds',
        vmin=2,
        vmax=6,
        colorbar_loc='right',
    )

    # Rasterize scatter dots (keeps axes/legend as vectors, dots as embedded raster)
    for coll in ax.collections:
        coll.set_rasterized(True)

    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('CEACAM6')

    # Remove frame
    for spine in ax.spines.values():
        spine.set_visible(False)

    cbar = ax.collections[0].colorbar
    if cbar:
        cbar.set_label('Expression')

    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    style.save_panel(fig, Path(OUTPUT_DIR) / "ceacam6_expression_umap")
    print(f"Saved: {Path(OUTPUT_DIR) / 'ceacam6_expression_umap'}"
          f".[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
