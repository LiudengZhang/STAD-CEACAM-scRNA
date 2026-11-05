#!/usr/bin/env python3
"""Shared scanpy dot plot for the supplementary marker panels.

  frames(adata, genes, groupby, order, path)  per-group scaled means
      (standard_scale='var') and fractions as a long table
      (group, gene, mean_scaled, fraction); raises if a gene is missing.
  draw(frames, order, genes, ...)  rebuilds the dot plot from that table via
      scanpy's dot_color_df / dot_size_df.

Run: python _dotplot.py --self-test
"""

from pathlib import Path

import numpy as np
import pandas as pd

if __name__ == "__main__":
    import sys as _s
    from pathlib import Path as _P
    _s.path[:0] = [str(_P(__file__).resolve().parents[3] / "00_Config"),
                   str(_P(__file__).resolve().parent)]

import panel_style_cns as style
from cnsfig import legend as cnslegend

def frames(adata, genes, groupby, order, *, path):
    """Per-group dot-plot values. `path` is the h5ad `adata` was read from; it
    selects the expression matrix (shared/expression.py)."""
    import scanpy as sc
    from shared.expression import (assert_x_stands_in_for_raw,
                                   expression_source)
    assert_x_stands_in_for_raw(adata, path)
    src = expression_source(adata, path)
    absent = [g for g in genes if g not in set(src.var_names)]
    if absent:
        raise SystemExit(
            f"_dotplot.frames: {', '.join(absent)} is not in the expression "
            f"matrix of {Path(path).name} ({len(src.var_names):,} genes). "
            f"Refusing to drop a declared gene: a marker plot missing a "
            f"marker looks exactly like one that never declared it. "
            f"Either the symbol is wrong, or the object is, or the gene has "
            f"to come out of the script's MARKERS.")
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
         cbar_title="Mean expression\nin group", panel="", shared_key=None,
         var_groups=None):
    """Draw the dot plot; returns (fig, main_ax).

    `var_groups`: ((start, end), label) pairs for scanpy's gene-group brackets.
    `shared_key`: omit the key (another panel's key decodes this one).
    """
    import anndata as ad
    import scanpy as sc
    missing_g = [g for g in genes if g not in set(fr["gene"])]
    missing_o = [o for o in order if o not in set(fr[groupby])]
    if missing_g or missing_o:
        raise SystemExit(
            "_dotplot.draw: the cached table does not cover what this panel "
            "declares"
            + (f"\n  genes missing: {', '.join(missing_g)}" if missing_g else "")
            + (f"\n  groups missing: {', '.join(missing_o)}" if missing_o else "")
            + "\nRefusing to draw a smaller plot than the script asks for. "
              "Run the panel with --recompute; if the gene is genuinely not in "
              "the matrix, frames() will now say so by name.")
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
                            dot_color_df=color_df, dot_size_df=size_df,
                            **({} if not var_groups else dict(
                                var_group_positions=[g[0] for g in var_groups],
                                var_group_labels=[g[1] for g in var_groups])))
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
    # Place the matrix by hand between the left margin and the key column.
    pos = main_ax.get_position()
    main_ax.set_position([left_mm / w_mm, pos.y0,
                          (w_mm - left_mm - key_mm - key_gap_mm) / w_mm, pos.height])
    if var_groups:
        # The bracket axes is scanpy's, above the matrix; it has to follow the
        # matrix across or the brackets stand over the wrong columns.
        gax = axes.get('gene_group_ax')
        if gax is None:
            raise RuntimeError(f"{panel}: var_groups given but scanpy drew no "
                               f"gene_group_ax")
        mp, gp = main_ax.get_position(), gax.get_position()
        gax.set_position([mp.x0, gp.y0, mp.width, gp.height])
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
    return fig, main_ax


# ---------------------------------------------------------------------------
# Self-test: draw() must refuse a cached table missing a gene or a group.
# ---------------------------------------------------------------------------

def _self_test():
    import pandas as _pd

    order = ["A", "B"]
    genes = ["G1", "G2", "G3"]
    full = _pd.DataFrame([{"group": g, "gene": v, "mean_scaled": 0.5,
                           "fraction": 0.5}
                          for g in order for v in genes])
    fails = []

    def case(label, fr, want):
        try:
            draw(fr, order, genes, w_mm=40, h_mm=20, left_mm=10, bottom_mm=10)
            got = None
        except SystemExit as exc:
            got = str(exc)
        except Exception as exc:                            # noqa: BLE001
            got = f"__unexpected__ {type(exc).__name__}: {exc}"
        if want is None:
            ok = got is None or not got.startswith("_dotplot.draw")
            print(f"  {label:44s} {'ok' if ok else 'FAIL'}")
            if not ok:
                fails.append(f"{label}: refused a complete frame: {got}")
            return
        ok = got is not None and want in got
        print(f"  {label:44s} "
              f"{'refused: yes' if ok else 'NOT CAUGHT / WRONG REASON'}")
        if not ok:
            fails.append(f"{label}: {got}")

    print("_dotplot self-test")
    case("a complete frame draws", full, None)
    case("MUTANT one declared gene missing from the cache",
         full[full["gene"] != "G2"], "genes missing: G2")
    case("MUTANT one declared group missing from the cache",
         full[full["group"] != "B"], "groups missing: B")

    if fails:
        for f in fails:
            print("   ", f)
        return 1
    print("  all cases pass")
    return 0


if __name__ == "__main__":
    import sys as _sys
    if "--self-test" in _sys.argv:
        raise SystemExit(_self_test())
    raise SystemExit("_dotplot.py is a driver; run a panel script, or "
                     "--self-test")
