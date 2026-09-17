#!/usr/bin/env python3
"""
One scanpy dot plot, drawn at print size from a cached table.

The four supplementary dot plots (S2 A, S2 B, S3 A, S4 A) were each a
`sc.pl.dotplot(adata, ...)` call on a 4x canvas that the page then fitted at
0.2-0.3, and each one opened a 12 GB h5ad to draw. This module splits them
the way Figure 4 E and 5 H were split on 2026-09-15:

  frames(adata, genes, groupby, order)   the computing half - scanpy's own
                                         DotPlot summarises the means
                                         (standard_scale='var') and the
                                         fractions, returned as a long table
                                         (group, gene, mean_scaled, fraction)
  draw(frames, order, genes, ...)        the drawing half - scanpy's DotPlot
                                         built from the two frames through its
                                         documented dot_color_df / dot_size_df
                                         inputs on a scaffold AnnData that
                                         carries the category order and no
                                         expression

The colour map, the scaling and the dot sizes RELATIVE to one another are
scanpy's, as before; the largest dot is set in millimetres so that dots at a
2.6 mm gene pitch do not run into one another, and every dot scales by the
same factor (compare_panel_content reads sizes as ratios). The size key and
the colour bar are drawn by cnsfig.legend.scanpy_compact_key in a reserved
column, as Figure 4 E's are.

ONE KEY FOR TWO DOT PLOTS ON ONE ROW (`shared_key`), since the evening of
2026-09-15. S2 A and S2 B share a row so that S1's second page fits the
main-figure page height; every dot plot drawn here is scaled alike - means
standard-scaled per gene to 0-1 on the same colour map, fractions on the
same largest dot - so one key decodes both, and F is drawn without a key
column, naming the panel whose key it reads (`shared_key="S2 B"`). That is
the only way round cnsfig.legend.require_size_key, and it is written into
the drawing's own call, not into the check.
"""

import numpy as np
import pandas as pd

import panel_style_cns as style
from cnsfig import legend as cnslegend


def frames(adata, genes, groupby, order):
    """The computing half. `genes` not in adata.var_names are dropped, in
    order - the predecessors filtered on the processed matrix's genes and
    then read expression from .raw, and so does this."""
    import scanpy as sc
    genes = [g for g in genes if g in adata.var_names]
    dp = sc.pl.DotPlot(adata, var_names=genes, groupby=groupby,
                       categories_order=order,
                       use_raw=(adata.raw is not None), standard_scale='var')
    color, size = dp.dot_color_df, dp.dot_size_df
    rows = [{"group": g, "gene": v, "mean_scaled": float(color.loc[g, v]),
             "fraction": float(size.loc[g, v])} for g in order for v in genes]
    return pd.DataFrame(rows)


def draw(fr, order, genes, *, w_mm, h_mm, left_mm, bottom_mm, top_mm=2.0,
         key_mm=16.5, key_gap_mm=2.5, largest_mm=2.2, rotation=90,
         italic=True, groupby="group",
         size_title="Fraction of cells\nin group (%)",
         cbar_title="Mean expression\nin group", panel="", shared_key=None):
    """The drawing half. Returns (fig, main_ax). `shared_key` names the
    panel on the same row whose key decodes this plot; the key column is
    then not reserved and no key is drawn here (see the module docstring)."""
    import anndata as ad
    import scanpy as sc
    genes = [g for g in genes if g in set(fr["gene"])]
    color_df = fr.pivot(index="group", columns="gene", values="mean_scaled") \
                 .loc[order, genes]
    size_df = fr.pivot(index="group", columns="gene", values="fraction") \
                .loc[order, genes]
    scaffold = ad.AnnData(
        X=np.zeros((len(order), len(genes)), dtype=np.float32),
        obs=pd.DataFrame({groupby: pd.Categorical(order, categories=order)},
                         index=[f"scaffold_{i}" for i in range(len(order))]),
        var=pd.DataFrame(index=genes))
    dotplot = sc.pl.dotplot(scaffold, var_names=genes, groupby=groupby,
                            categories_order=order, use_raw=False, cmap='Reds',
                            show=False, save=None, return_fig=True,
                            figsize=style.figsize_mm(w_mm, h_mm),
                            dot_color_df=color_df, dot_size_df=size_df)
    # scatter's s is the diameter squared in points; cmap is passed again
    # because style() rewrites every field it takes.
    dotplot.style(largest_dot=(largest_mm * style.PT_PER_MM) ** 2, cmap='Reds')
    dotplot.legend(width=key_mm / 25.4)
    dotplot.make_figure()
    axes = dotplot.get_axes()
    main_ax = axes['mainplot_ax']
    main_ax.set_xticklabels([t.get_text() for t in main_ax.get_xticklabels()],
                            rotation=rotation, ha='right' if rotation < 90 else 'center',
                            style='italic' if italic else 'normal',
                            fontsize=style.tick_pt())
    main_ax.tick_params(labelsize=style.tick_pt())
    for ax in dotplot.fig.get_axes():
        for coll in ax.collections:
            coll.set_rasterized(False)
        for img in ax.images:
            img.set_rasterized(False)
    fig = dotplot.fig
    if shared_key:
        key_mm, key_gap_mm = 0.0, 0.0
    style.margins_mm(fig, left=left_mm, right=key_mm + key_gap_mm, top=top_mm,
                     bottom=bottom_mm)
    # scanpy's gridspec keeps a legend column of its own; the matrix is placed
    # by hand from the left margin to key_gap_mm short of the key column.
    pos = main_ax.get_position()
    main_ax.set_position([left_mm / w_mm, pos.y0,
                          (w_mm - left_mm - key_mm - key_gap_mm) / w_mm, pos.height])
    if shared_key:
        for name in ("size_legend_ax", "color_legend_ax"):
            ax = axes.get(name)
            if ax is not None:
                ax.remove()
        print(f"  {panel}: no key of its own - decoded by the key of {shared_key} "
              f"on the same row (same standard_scale, same largest dot)")
    else:
        key, areas = cnslegend.scanpy_compact_key(
            dotplot, axes, key_mm, size_title=size_title, cbar_title=cbar_title,
            label_pt=style.tick_pt(), title_pt=style.tick_pt())
        cnslegend.require_size_key(key, dot_areas=areas, panel=panel)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(f"{panel}: ink outside the {w_mm} x {h_mm} mm canvas "
                           f"(l,r,b,t mm): {over}")
    return fig, main_ax
