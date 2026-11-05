#!/usr/bin/env python3
"""UMAP coloured by a per-cell value, drawn from a cached (x, y, value) table.
Shared by Supplementary panels S5 B, C and D through `panel()`; cells are
plotted in ascending order of value, with a colour bar beside the map.
"""

from pathlib import Path

import numpy as np
import pandas as pd

import panel_style_cns as style

FIG = "S5_CD8_TCells"
#: Panel geometry in mm.
W_MM, H_MM, MAP_MM, LEFT_MM, TOP_MM = 53.0, 46.0, 36.0, 4.5, 3.6

def frame(adata, values):
    """(x, y, value) for every cell of `adata`, in its order."""
    xy = np.asarray(adata.obsm["X_umap"], dtype=float)
    return pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1],
                         "value": np.asarray(values, dtype=float)})

def draw(fr, *, title, cmap, vmin=None, vmax=None, w_mm, h_mm, map_mm,
         left_mm, top_mm, cbar_w_mm=1.5, cbar_gap_mm=1.5, dot_mm=0.25):
    """Returns (fig, ax). Axes are created map first, colour bar second."""
    fig = style.figure_mm(w_mm, h_mm)
    ax = fig.add_axes([left_mm / w_mm, 1 - (top_mm + map_mm) / h_mm,
                       map_mm / w_mm, map_mm / h_mm])
    v = fr["value"].to_numpy()
    order = np.argsort(-v, kind="stable")[::-1]
    sc_ = ax.scatter(fr["x"].to_numpy()[order], fr["y"].to_numpy()[order], c=v[order],
                     cmap=cmap, vmin=vmin, vmax=vmax, s=(dot_mm * style.PT_PER_MM) ** 2,
                     edgecolors="none", linewidths=0, rasterized=True)
    ax.set_title(title, pad=2)
    ax.set_xlabel("UMAP1", fontsize=style.tick_pt())
    ax.set_ylabel("UMAP2", fontsize=style.tick_pt())
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(style.RULE_PT)
    cax = fig.add_axes([(left_mm + map_mm + cbar_gap_mm) / w_mm,
                        1 - (top_mm + 0.9 * map_mm) / h_mm, cbar_w_mm / w_mm,
                        0.8 * map_mm / h_mm])
    cb = fig.colorbar(sc_, cax=cax)
    cb.ax.tick_params(labelsize=style.tick_pt(), width=style.RULE_PT, length=1.5)
    cb.outline.set_linewidth(style.EDGE_PT)
    return fig, ax

def stomach(adata):
    """The stomach cells."""
    return adata[adata.obs['Sample site'] == 'Stomach'].copy()

def gene_values(adata, gene, *, path):
    """One gene's expression per cell. `path` is the h5ad `adata` was read
    from; it selects the expression matrix (shared/expression.py)."""
    from shared.expression import (assert_x_stands_in_for_raw,
                                   expression_source)
    assert_x_stands_in_for_raw(adata, path)
    src_ = expression_source(adata, path)
    col = src_[:, gene].X
    return col.toarray().ravel() if hasattr(col, 'toarray') else col.ravel()

def panel(here, panel_dir, compute, *, title, cmap, vmin=None):
    """Draw and save one S5 UMAP panel; `compute` runs only if data/cells.csv
    is missing or with --recompute."""
    import _driver_base as base
    from cnsfig import cache
    base.apply_style()
    fr = cache.table(Path(here), "cells", compute)
    fig, _ = draw(fr, title=title, cmap=cmap, vmin=vmin, w_mm=W_MM, h_mm=H_MM,
                  map_mm=MAP_MM, left_mm=LEFT_MM, top_mm=TOP_MM)
    base.save(fig, FIG, panel_dir, f"panel_{panel_dir}")
    return 0
