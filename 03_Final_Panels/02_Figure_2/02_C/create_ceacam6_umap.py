#!/usr/bin/env python3
"""
Figure 2 panel C - CEACAM6 expression on the epithelial UMAP.

The panel is drawn at the millimetre rectangle it prints in, read from
03_Final_Panels/panel_rects.csv through 00_Config/slots.py, so the type size
set here is the type size printed. The margins are millimetres of paper and are
set before the colour bar is made, because scanpy sizes the bar from the axes
box it finds; they are not fitted to the ink afterwards, since moving the
subplot parameters after the bar exists would leave the bar behind.

  printed panel  Figure 2 C       (PROVENANCE.csv; the directory name agrees,
                                   but the letter was looked up, not inferred)

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the colour bar's tick labels - at 5 * SCALE. MARK
    carries the non-type point sizes across to the 1:1 canvas:

        MARK = style.tick_pt() / (SMALL_PT * SCALE)
        AREA = MARK ** 2

    The point area takes AREA. The spine and tick widths the earlier drawing
    set are not carried over: those are axes furniture, cnsplots has its own
    settings for them, and following the library rather than rescaling the old
    numbers is the standard-methods rule.

    The top margin keeps the axes below the corner the panel letter is printed
    in: the points cover their whole axes box, so the corner can only be left
    clear by placing the box below it.

Every cell, every colour scale limit and every string is the earlier drawing's.
The drawing code is the same code.
"""

import scanpy as sc
import numpy as np
import os

# Central config
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *                                       # noqa: E402,F403
import panel_style_cns as style                           # noqa: E402
import slots                                              # noqa: E402

SCALE = 4                           # the earlier canvas multiplier
SMALL_PT = 5.0                      # the earlier smallest body type

PANEL_LETTER = "C"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(2, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(2, PANEL_LETTER)
MARGIN = dict(left=2.6, right=7.0, top=4.2, bottom=1.5)

# Paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = BASE_DIR


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
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
    # The gene symbol is set in italic, as the shipped page sets it and as
    # the rest of this figure sets one. Style only; the string is unchanged.
    ax.set_title('CEACAM6', fontstyle='italic')

    # Remove frame
    for spine in ax.spines.values():
        spine.set_visible(False)

    cbar = ax.collections[0].colorbar
    if cbar:
        cbar.set_label('Expression')

    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    style.save_panel(fig, Path(OUTPUT_DIR) / "ceacam6_expression_umap")
    print(f"Saved: {Path(OUTPUT_DIR) / 'ceacam6_expression_umap'}"
          f".[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
