#!/usr/bin/env python3
"""
Figure 2B: CEACAM5 expression on the epithelial UMAP. draw() is shared with
Figure 2C (02_C/create_ceacam6_umap.py), which calls it for CEACAM6.

Reads the epithelial h5ad (EPITHELIAL_H5AD); log-normalised expression from
.raw, or from .X when the file has no .raw.
Writes ceacam5_expression_umap.[svg|pdf|png] beside this script.

Run: python create_ceacam5_umap.py

The colour scale is fixed at 2-6 for both genes. Margins are set before scanpy
adds the colour bar, because the bar is sized from the axes box.
"""

import scanpy as sc
import numpy as np

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import EPITHELIAL_H5AD                           # noqa: E402
from shared.expression import assert_x_stands_in_for_raw    # noqa: E402
from cnsfig import Panel, style                             # noqa: E402

SCALE = 4                           # sets the point size relative to the tick size
SMALL_PT = 5.0

MARGIN = dict(left=2.6, right=7.0, top=4.2, bottom=1.5)

# Each calling script builds its own Panel, so the output lands in its directory.
PANEL = Panel(2, "B", __file__, small_pt=SMALL_PT)


def draw(panel, gene, stem):
    """One gene's expression UMAP, at the slot `panel`'s printed letter prints in."""
    style.apply(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2

    adata = sc.read_h5ad(EPITHELIAL_H5AD)
    # Raises if this file has no .raw and its .X is not known to equal it.
    assert_x_stands_in_for_raw(adata, EPITHELIAL_H5AD)

    # The deposited h5ad has no .raw; its .X holds the same log1p matrix.
    src = adata.raw if adata.raw is not None else adata
    idx = src.var_names.get_loc(gene)
    expr = src.X[:, idx]
    if hasattr(expr, 'toarray'):
        expr = expr.toarray().flatten()
    else:
        expr = np.array(expr).flatten()
    color = f'{gene}_expr'
    adata.obs[color] = expr

    fig, ax = style.subplots_mm(panel.w_mm, panel.h_mm)
    # Before the colour bar (see header).
    style.margins_mm(fig, **MARGIN)

    sc.pl.umap(
        adata,
        color=color,
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

    # Dots as an embedded raster; the type and the bar stay vectors.
    for coll in ax.collections:
        coll.set_rasterized(True)

    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks([])
    ax.set_yticks([])
    # Gene symbol in italic.
    ax.set_title(gene, fontstyle='italic')

    for spine in ax.spines.values():
        spine.set_visible(False)

    cbar = ax.collections[0].colorbar
    if cbar:
        cbar.set_label('Expression')

    panel.finish(fig, stem)


def main():
    draw(PANEL, 'CEACAM5', "ceacam5_expression_umap")


if __name__ == '__main__':
    main()
