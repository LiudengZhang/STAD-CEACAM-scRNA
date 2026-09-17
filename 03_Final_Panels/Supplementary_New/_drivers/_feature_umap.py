#!/usr/bin/env python3
"""
One UMAP coloured by a per-cell value, drawn at print size from a cached table.

Supplementary panels S4 B, C and D were `sc.pl.umap(adata, color=...)` on a
4x canvas. The drawing half here reproduces what scanpy 1.9.6 draws - the
cells in ascending order of the value (`np.argsort(-v, kind="stable")[::-1]`),
a square scatter with the value mapped through `cmap` between `vmin` and
`vmax`, a colour bar beside it - from a (x, y, value) table, so the h5ad is
not reopened to redraw and compare_panel_content.py sees the same offsets in
the same order.
"""

import numpy as np
import pandas as pd

import panel_style_cns as style


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
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(f"{title}: ink outside the {w_mm} x {h_mm} mm canvas: {over}")
    return fig, ax
